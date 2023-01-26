/// @file tandem_dc_motors.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <exception>
#include <iostream>
#include <iomanip>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "solver/detail/observer.hpp"
#include "solver/stepper/stepper_adaptive.hpp"
#include "solver/stepper/stepper_euler.hpp"
#include "solver/stepper/stepper_rk4.hpp"
#include "solver/solver.hpp"

#include "defines.hpp"

namespace tandem_dc_motors
{

#define AMBIENT_TEMPERATURE 22.

/// @brief State of the system.
///     x[0] : Current motor 1.
///     x[1] : Current motor 2.
///     x[2] : Angular speed motor 1.
///     x[3] : Angular speed motor 2.
///     x[4] : Angular speed load.
///     x[5] : Angle shaft 1.
///     x[6] : Angle shaft 2.
///     x[7] : The temperature of motor 1.
///     x[8] : The temperature of motor 2.
using State = std::array<Variable, 9>;

struct Parameters {
    /// Supplied voltage motor 1 [V].
    Variable v_a1;
    /// Winding resistance motor 1 [Ohms].
    Variable R_a1;
    /// Winding inductance motor 1 [Henrys].
    Variable L_a1;
    /// Back-EMF motor 1 [V * s / rad].
    Variable K_e1;
    /// Torque constant motor 1 [N * m / A].
    Variable Kt_m1;
    /// Angular momentum motor 1 [kg.m ^ 2].
    Variable J_m1;
    /// Coulomb friction motor 1 [N.m].
    Variable Kd_m1;
    /// Thermal resistance of motor 1 [C / Watt].
    Variable R_Th_m1;
    /// Thermal capacity of the coil of motor 1 [Joule / C].
    Variable C_Th_m1;

    /// Supplied voltage motor 2 [V].
    Variable v_a2;
    /// Winding resistance motor 2 [Ohms].
    Variable R_a2;
    /// Winding inductance motor 2 [Henrys].
    Variable L_a2;
    /// Back-EMF motor 2 [V * s / rad].
    Variable K_e2;
    /// Torque constant motor 2 [N * m / A].
    Variable Kt_m2;
    /// Angular momentum motor 2 [kg.m ^ 2].
    Variable J_m2;
    /// Coulomb friction motor 2 [N.m].
    Variable Kd_m2;
    /// Thermal resistance of motor 2 [C / Watt].
    Variable R_Th_m2;
    /// Thermal capacity of the coil of motor 2 [Joule / C].
    Variable C_Th_m2;

    /// Coulombic friction of shaft 1.
    Variable Kd_s1;
    /// Elasticity of shaft 1.
    Variable Ke_s1;

    /// Coulombic friction of shaft 2.
    Variable Kd_s2;
    /// Elasticity of shaft 2.
    Variable Ke_s2;

    /// Load torque.
    Variable T_l;
    /// Load inertia.
    Variable J_l;
    /// Coulombic friction of load.
    Variable Kd_l;

    // Motor 1/Motor 2 gear ratio.
    Variable Gr;
    /// Ambient temperature.
    Variable T_Amb;

    Parameters()
        : // Motor 1
          v_a1(12),
          R_a1(12.00),
          L_a1(18.00e-03),
          K_e1(21.22),
          Kt_m1(0.059),
          J_m1(0.5),
          Kd_m1(0.02),
          R_Th_m1(2.2),
          C_Th_m1(9 / R_Th_m1),
          // Motor 2
          v_a2(12),
          R_a2(12.00),
          L_a2(18.00e-03),
          K_e2(21.22),
          Kt_m2(0.059),
          J_m2(0.5),
          Kd_m2(0.02),
          R_Th_m2(2.2),
          C_Th_m2(9 / R_Th_m2),
          // Shaft 1
          Kd_s1(0.05),
          Ke_s1(0.01),
          // Shaft 2
          Kd_s2(0.05),
          Ke_s2(0.01),
          // Load
          T_l(this->compute_load_torque(0.2)),
          J_l(this->compute_inertia_load(0.2)),
          Kd_l(0.05),
          // Gear ratio
          Gr(1 / 2),
          // Ambient temperature
          T_Amb(AMBIENT_TEMPERATURE)
    {
        // Nothing to do.
    }

private:
    /// @brief Copmutes the torque for the load.
    /// @param m_l Load mass [kg].
    /// @param F_A External Force [N].
    /// @param r Radius of the roller [m].
    /// @param mu Friction coefficient.
    /// @param theta Inclination angle [rad].
    /// @param eta Efficiency.
    /// @param G_r Gear ratio.
    /// @return the torque for the load.
    inline Variable compute_load_torque(
        const Variable m_l   = 1,
        const Variable F_A   = 0,
        const Variable r     = 0.20,
        const Variable mu    = 1,
        const Variable theta = 0,
        const Variable eta   = 1,
        const Variable G_r   = 1)
    {
        // Gravitational force [N].
        constexpr Variable g = 9.8;
        // Diameter of the roller.
        const Variable D = r * 2;

        Variable F = F_A + m_l * g * (std::sin(theta) + mu * std::cos(theta));
        return (F * D) / (2 * eta * G_r);
    }

    /// @brief Compute the inertial of the load (object + conveyor belt).
    /// @param m_l Mass of the load [kg].
    /// @param m_r Mass of the roller [kg].
    /// @param m_b Mass of the belt [kg].
    /// @param r Radius of the roller [m].
    /// @return Variable
    inline Variable compute_inertia_load(
        const Variable m_l = 1,
        const Variable m_r = 0.5,
        const Variable m_b = 1,
        const Variable r   = 0.20)
    {
        /// JA : Moment of inertia of the roller [kg・m2]
        const Variable J_r = this->compute_inertia_cylinder(m_r, r);
        // Diameter of the roller.
        const Variable D = r * 2;
        // Compute the inertial of the conveyor belt.
        return J_r + ((m_l + m_b) * D * D) / 4;
    }

    /// @brief Compute the moment of inertia of the central axis of a cylinder.
    /// @param m Mass of cylinder [kg].
    /// @param r Radius of the cylinder [m].
    /// @return the moment of inertia.
    inline Variable compute_inertia_cylinder(const Variable m, const Variable r)
    {
        return (m * r * r) / 2;
    }
};

/// @brief The model.
struct Model : public Parameters {
    Model(Parameters _param = Parameters())
        : Parameters(std::move(_param))
    {
        // Nothing to do.
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    constexpr inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        const auto i_m1 = x[0]; // Current motor 1.
        const auto i_m2 = x[1]; // Current motor 2.
        const auto w_m1 = x[2]; // Angular speed motor 1.
        const auto w_m2 = x[3]; // Angular speed motor 2.
        const auto w_l  = x[4]; // Angular speed load.
        const auto a_s1 = x[5]; // Angle shaft 1.
        const auto a_s2 = x[6]; // Angle shaft 2.
        const auto t_m1 = x[7]; // The temperature of motor 1.
        const auto t_m2 = x[8]; // The temperature of motor 2.

        // Current motor 1.
        dxdt[0] = v_a1                // Voltage source
                - R_a1 / L_a1 * i_m1  // Resistance
                - K_e1 / L_a1 * w_m1; // Back EMF

        // Current motor 2.
        dxdt[1] = v_a2                // Voltage source
                - R_a2 / L_a2 * i_m2  // Resistance
                - K_e2 / L_a2 * w_m2; // Back EMF

        // Angular speed motor 1.
        dxdt[2] = Kt_m1 / J_m1 * i_m1                // Torque applied to M1
                - Gr * (Kd_m1 + Kd_s1) / J_m1 * w_m1 // Impact of static friction from M1 itself, and S1 based on the speed of M1
                + Gr * Kd_s1 / J_m1 * w_m2           // Impact of the coulombic friction of S1 based on the speed of M2
                - Gr * Ke_s1 / J_m1 * a_s1;          // Impact of the elasticity of S1 based on the angle of S1

        // Angular speed motor 2.
        dxdt[3] = Kt_m2 / J_m2 * i_m2                        // Torque applied to M2
                + Gr * Kd_s1 / J_m2 * w_m1                   // Impact of the coulombic friction of S2 based on the speed of M1
                - Gr * (Kd_m1 + Kd_m2 + Kd_s1) / J_m2 * w_m2 // Impact of static friction from M2 itself, M1, and S1 based on the speed of M2
                + Gr * Kd_s2 / J_m2 * w_l                    // Impact of the coulombic friction of S2 based on the speed of the load
                + Ke_s1 / J_m2 * a_s1                        // Impact of the elasticity of S1 based on the angle of S1
                - Ke_s2 / J_m2 * a_s2;                       // Impact of the elasticity of S2 based on the angle of S2

        // Angular speed load.
        dxdt[4] = Gr * Kd_s2 / J_l * w_m2         // Impact of the coulombic friction of S2 based on the speed of M2
                - Gr * (Kd_l + Kd_s2) / J_l * w_l // Impact of static friction from the load itself, and S2 based on the speed of the load
                + Ke_s2 / J_l * a_s2              // Impact of the elasticity of S2 based on the angle of S2
                + T_l;                            // The torque of the load

        // Angle shaft 1.
        dxdt[5] = w_m1 - w_m2; // The angle of the first shaft computed as the difference in speed between M1 and M2

        // Angle shaft 2.
        dxdt[6] = w_m2 - w_l; // The angle of the first shaft computed as the difference in speed between M2 and the load

        // The temperature of motor 1.
        dxdt[7] = (R_a1 / C_Th_m1) * i_m1 * i_m1 + (T_Amb - t_m1) / (C_Th_m1 * R_Th_m1);

        // The temperature of motor 2.
        dxdt[8] = (R_a2 / C_Th_m2) * i_m2 * i_m2 + (T_Amb - t_m2) / (C_Th_m2 * R_Th_m2);
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, i_a1, i_a2, w_m1, w_m2, w_l, a_s1, a_s2, t_m1, t_m2;

    ObserverSave() = default;

    constexpr inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            i_a1.emplace_back(x[0]);
            i_a2.emplace_back(x[1]);
            w_m1.emplace_back(x[2]);
            w_m2.emplace_back(x[3]);
            w_l.emplace_back(x[4]);
            a_s1.emplace_back(x[5]);
            a_s2.emplace_back(x[6]);
            t_m1.emplace_back(x[7]);
            t_m2.emplace_back(x[8]);
        }
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverPrint : public solver::detail::DecimationObserver<DECIMATION> {
    ObserverPrint() = default;
    inline void operator()(const State &x, const Time &t)
    {
        if (this->observe())
            std::cout << std::fixed << std::setprecision(4) << t << " " << x << "\n";
    }
};

} // namespace tandem_dc_motors

int main(int, char **)
{
    using namespace tandem_dc_motors;

    // Instantiate the model.
    Model model;

    // Initial and runtime states.
    const State x0{ .0, .0, .0, .0, .0, .0, .0, AMBIENT_TEMPERATURE, AMBIENT_TEMPERATURE };
    State x;

    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 100.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = solver::ErrorFormula::Mixed;
    const auto Iterations = 3;
    using Rk4             = solver::stepper_rk4<State, Time>;
    using AdaptiveRk4     = solver::stepper_adaptive<Rk4, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveRk4 adaptive_rk4(time_delta);

    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    ObserverSave<0> obs_adaptive_rk4;
#elif 1
    ObserverPrint<0> obs_adaptive_rk4;
#endif

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";
    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    solver::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw.partials()[0] << "\n";

#ifdef SC_ENABLE_PLOT
    matplot::hold(matplot::on);
    // matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.i_a1, "--g")->line_width(2).display_name("Current M1");
    // matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.i_a2, "-.g")->line_width(2).display_name("Current M2");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_m1, "--b")->line_width(2).display_name("Angular Speed M1");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_m2, "-.b")->line_width(2).display_name("Angular Speed M2");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.a_s1, "--r")->line_width(2).display_name("Angle Shaft 1");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.a_s2, "-.r")->line_width(2).display_name("Angle Shaft 2");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.t_m1, "--m")->line_width(2).display_name("Temperature Motor 1");
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.t_m2, "-.m")->line_width(2).display_name("Temperature Motor 2");
    // matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_l)->line_width(2).display_name("Angular Speed Load");
    matplot::legend(matplot::on);
    matplot::show();
#endif
    return 0;
}