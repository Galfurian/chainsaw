/// @file tandem_dc_motors.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <exception>
#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

namespace tandem_dc_motors
{

/// Ambient temperature.
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

/// Defines the mode of operation of the model.
enum Mode {
    mode_0, ///< Va1 = 12V, Va2 = 12V
    mode_1, ///< Va1 = 24V, Va2 = 24V
    mode_2, ///< Va1 = 36V, Va2 = 36V
    mode_3, ///< Va1 = 48V, Va2 = 48V
    mode_4, ///< Va1 = 60V, Va2 = 60V
};

/// Controls the mode of operation of the model.
struct Step {
    Mode mode;     ///< Which mode we need to use in a given sequence.
    Time duration; ///< How long the sequence will last, in seconds.
};

/// Defines the sequence of modes, we need to use during the simulation.
using Sequence = std::vector<Step>;

/// Parameters of our model.
struct Parameters {
    /// Current operating mode.
    Mode mode;
    /// Supplied voltage motor 1 [V].
    Variable v_a1;
    /// Winding resistance motor 1 [Ohms].
    Variable Ra_m1;
    /// Winding inductance motor 1 [Henrys].
    Variable La_m1;
    /// Back-EMF motor 1 [V * s / rad].
    Variable Ke_m1;
    /// Torque constant motor 1 [N * m / A].
    Variable Kt_m1;
    /// Angular momentum motor 1 [kg.m ^ 2].
    Variable J_m1;
    /// Coulomb friction motor 1 [N.m].
    Variable Kd_m1;
    /// Thermal resistance of motor 1 [C / Watt].
    Variable Rth_m1;
    /// Thermal capacity of the coil of motor 1 [Joule / C].
    Variable Cth_m1;

    /// Supplied voltage motor 2 [V].
    Variable v_a2;
    /// Winding resistance motor 2 [Ohms].
    Variable Ra_m2;
    /// Winding inductance motor 2 [Henrys].
    Variable La_m2;
    /// Back-EMF motor 2 [V * s / rad].
    Variable Ke_m2;
    /// Torque constant motor 2 [N * m / A].
    Variable Kt_m2;
    /// Angular momentum motor 2 [kg.m ^ 2].
    Variable J_m2;
    /// Coulomb friction motor 2 [N.m].
    Variable Kd_m2;
    /// Thermal resistance of motor 2 [C / Watt].
    Variable Rth_m2;
    /// Thermal capacity of the coil of motor 2 [Joule / C].
    Variable Cth_m2;

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

    /// Ambient temperature.
    Variable T_Amb;
    /// Temperature coefficient of motor winding (1/K).
    Variable tco_wind_m1;
    /// Temperature coefficient of permanent magnet (1/K).
    Variable tco_magn_m1;
    /// Temperature coefficient of motor winding (1/K).
    Variable tco_wind_m2;
    /// Temperature coefficient of permanent magnet (1/K).
    Variable tco_magn_m2;

    /// @brief Creates a new default parameter set.
    Parameters()
        : // Current mode.
          mode(mode_0),
          // Motor 1
          v_a1(12),
          Ra_m1(12),
          La_m1(800e-03),
          Ke_m1(21.22),
          Kt_m1(0.059),
          J_m1(0.5),
          Kd_m1(0.02),
          Rth_m1(2.2),
          Cth_m1(9 / Rth_m1),
          // Motor 2
          v_a2(12),
          Ra_m2(12),
          La_m2(800e-03),
          Ke_m2(21.22),
          Kt_m2(0.059),
          J_m2(0.5),
          Kd_m2(0.02),
          Rth_m2(2.2),
          Cth_m2(9 / Rth_m2),
          // Shaft 1
          Kd_s1(0.05),
          Ke_s1(0.01),
          // Shaft 2
          Kd_s2(0.05),
          Ke_s2(0.01),
          // Load
          T_l(this->compute_load_torque(0.2, 0, 0.10, 0.05)),
          J_l(this->compute_inertia_load(0.2, 0.2, 0.2, 0.10)),
          Kd_l(0.05),
          // Ambient temperature
          T_Amb(AMBIENT_TEMPERATURE),
          tco_wind_m1(0.004),
          tco_magn_m1(-0.002),
          tco_wind_m2(0.004),
          tco_magn_m2(-0.002)
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
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
    {
(void) t;
        // Change the input voltage, base on the mode.
        if (mode == mode_0)
            v_a1 = 12, v_a2 = 12;
        else if (mode == mode_1)
            v_a1 = 24, v_a2 = 24;
        else if (mode == mode_2)
            v_a1 = 36, v_a2 = 36;
        else if (mode == mode_3)
            v_a1 = 48, v_a2 = 48;
        else if (mode == mode_4)
            v_a1 = 60, v_a2 = 60;
        else
            std::exit(1);

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
        dxdt[0] = (v_a1            // Voltage source
                   - Ra_m1 * i_m1  // Resistance
                   - Ke_m1 * w_m1) // Back EMF
                / La_m1;

        // Current motor 2.
        dxdt[1] = (v_a2            // Voltage source
                   - Ra_m2 * i_m2  // Resistance
                   - Ke_m2 * w_m2) // Back EMF
                / La_m2;

        // Angular speed motor 1.
        dxdt[2] = (Kt_m1 * i_m1             // Torque applied to M1
                   - (Kd_m1 + Kd_s1) * w_m1 // Impact of static friction from M1 itself, and S1 based on the speed of M1
                   + Kd_s1 * w_m2           // Impact of the coulombic friction of S1 based on the speed of M2
                   - Ke_s1 * a_s1)          // Impact of the elasticity of S1 based on the angle of S1
                / J_m1;

        // Angular speed motor 2.
        dxdt[3] = (Kt_m2 * i_m2                     // Torque applied to M2
                   + Kd_s1 * w_m1                   // Impact of the coulombic friction of S2 based on the speed of M1
                   + Kd_s2 * w_l                    // Impact of the coulombic friction of S2 based on the speed of the load
                   - (Kd_m1 + Kd_m2 + Kd_s1) * w_m2 // Impact of static friction from M2 itself, M1, and S1 based on the speed of M2
                   + Ke_s1 * a_s1                   // Impact of the elasticity of S1 based on the angle of S1
                   - Ke_s2 * a_s2)                  // Impact of the elasticity of S2 based on the angle of S2
                / J_m2;

        // Angular speed load.
        dxdt[4] = (Kd_s2 * w_m2           // Impact of the coulombic friction of S2 based on the speed of M2
                   + Ke_s2 * a_s2         // Impact of the elasticity of S2 based on the angle of S2
                   - (Kd_l + Kd_s2) * w_l // Impact of static friction from the load itself, and S2 based on the speed of the load
                   + T_l)                 // The torque of the load
                / J_l;

        // Angle shaft 1.
        dxdt[5] = w_m1 - w_m2; // The angle of the first shaft computed as the difference in speed between M1 and M2

        // Angle shaft 2.
        dxdt[6] = w_m2 - w_l; // The angle of the first shaft computed as the difference in speed between M2 and the load

        // The temperature of motor 1.
        dxdt[7] = (Ra_m1 * i_m1 * i_m1        // Supplied Power.
                   - (t_m1 - T_Amb) / Rth_m1) // Power losses because of thermal dissipation.
                / Cth_m1;

        // The temperature of motor 2.
        dxdt[8] = (Ra_m2 * i_m2 * i_m2        // Supplied Power.
                   - (t_m2 - T_Amb) / Rth_m2) // Power losses because of thermal dissipation.
                / Cth_m2;

        Ra_m1 = Ra_m1 * (1 + tco_wind_m1 * (t_m1 - T_Amb));
        Kt_m1 = Kt_m1 * (1 + tco_magn_m1 * (t_m1 - T_Amb));

        Ra_m2 = Ra_m2 * (1 + tco_wind_m2 * (t_m2 - T_Amb));
        Kt_m2 = Kt_m2 * (1 + tco_magn_m2 * (t_m2 - T_Amb));
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
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
    std::vector<Variable> time, i_a1, i_a2, w_m1, w_m2, w_l, a_s1, a_s2, t_m1, t_m2;
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
    Time time             = 0.0;
    const Time time_start = 0.0;
    const Time time_end   = 300.0;
    const Time time_delta = 0.000001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    const auto Iterations = 10;
    using AdaptiveRk4     = chainsaw::stepper_adaptive<chainsaw::stepper_rk4<State, Time>, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveRk4 solver;
    solver.set_tollerance(1e-06);
    solver.set_min_delta(1e-09);
    solver.set_max_delta(1e-01);

    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs;

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";
    std::cout << "Simulating with `RK4`...\n";

    // Define the simulation sequence.
    Sequence sequence = {
        Step{ mode_0, 150. },
        Step{ mode_1, 150. },
        Step{ mode_2, 150. },
        Step{ mode_3, 150. },
        Step{ mode_4, 150. }
    };

    // Set the initial state.
    x = x0;
    // Start the simulation.
    sw.start();
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        // Set the mode.
        model.mode = sequence[i].mode;
        // Get the duration.
        Time duration = sequence[i].duration;
        // Run the solver.
        chainsaw::integrate_adaptive(solver, obs, model, x, time, time + duration, time_delta);
        // Advance time.
        time += duration;
    }
    // Get the elapsed time.
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive solver took " << std::setw(12) << solver.steps() << " steps, for a total of " << sw.partials()[0] << "\n";

#ifdef SC_ENABLE_PLOT
    auto figure = matplot::figure(true);
    // figure->position(0, 0, 800, 500);
    // figure->font_size(16);
    matplot::grid(matplot::on);
    matplot::hold(matplot::on);
    matplot::plot(obs.time, obs.w_m1, "-b")->line_width(1).display_name("Angular Speed M1 (rad/s)");
    matplot::plot(obs.time, obs.w_m2, "--b")->line_width(2).display_name("Angular Speed M2 (rad/s)");
    // matplot::plot(obs.time, obs.t_m1, "-m")->line_width(1).display_name("Temperature M1");
    // matplot::plot(obs.time, obs.t_m2, "--m")->line_width(2).display_name("Temperature M2");
    matplot::plot(obs.time, obs.i_a1, "-g")->line_width(1).display_name("Current M1 (A)");
    matplot::plot(obs.time, obs.i_a2, "--g")->line_width(2).display_name("Current M2 (A)");
    // matplot::plot(obs.time, obs.a_s1, "-r")->line_width(1).display_name("Angle Shaft 1");
    // matplot::plot(obs.time, obs.a_s2, "-.r")->line_width(2).display_name("Angle Shaft 2");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::xlabel("Time (s)");
    // matplot::ylabel("Temperature (C)");
    // matplot::show();

    char simulation_type[] = "speed";

    std::stringstream ss;
    ss << "result_dc_motor_" << simulation_type;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        ss << "_" << sequence[i].mode;
    }
    ss << ".png";
    matplot::save(ss.str());
#endif
    return 0;
}