/// @file dcmotor_v2.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-04-13

#include <timelib/stopwatch.hpp>
#include <exception>
#include <iostream>
#include <iomanip>

#ifdef ENABLE_PLOT
#include <gpcpp/gnuplot.hpp>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

namespace dcmotor_v2
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Temperature
using State = std::array<double, 3>;

/// @brief This one just containts the parameters.
struct Parameters {
    double E_0;    // [V]
    double Tau_L0; // [N.m]
    double T_Amb;  // [deg]
    double B_2C;   // [N]

    // motor parameters , Nachtigal , Table 16.5 p. 663
    double J_1; // in*oz*s ^2/ rad
    double B_1; // in*oz*s/ rad

    // electrical / mechanical relations
    double K_E; // back emf coefficient , e_m = K_E * omega_m ( K_E= alpha * omega )
    double K_T; // torque coeffic ., in English units K_T is not = K_E ! ( K_T = alpha * )
    double R_A; // Ohms
    double L_A; // H

    // gear - train and load parameters
    double J_2; // in*oz*s^2/ rad // 10x motor J
    double B_2; // in*oz*s/ rad ( viscous )
    double N;   // motor / load gear ratio ; omega_1 = N omega_2

    // Thermal model parameters
    double R_TM; // Thermal resistance (C / Watt)
    double C_TM; // Thermal capacity (Watt - sec / C) (-> 9 sec time constant - fast !)

    // Support variables.
    double Jeq;
    double Beq;
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;

    Parameters()
        : E_0(120.),
          Tau_L0(80.),
          T_Amb(18.),
          B_2C(300.),
          J_1(0.0035),
          B_1(0.064),
          K_E(0.1785),
          K_T(141.6 * K_E),
          R_A(8.4),
          L_A(0.0084),
          J_2(0.035),
          B_2(2.64),
          N(8.),
          R_TM(2.2),
          C_TM(9. / R_TM),
          Jeq(J_2 + N * 2 * J_1),
          Beq(B_2 + N * N * B_1),
          a(R_A / L_A),
          b(K_E * N / L_A),
          c(N * K_T / Jeq),
          d(Beq / Jeq),
          e(B_2C / Jeq),
          f(R_A / C_TM),
          g(1. / (C_TM * R_TM))
    {
        // Nothing to do.
    }
};

/// @brief The dc motor itself.
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
    constexpr void operator()(const State &x, State &dxdt, Time t) const noexcept
    {
        //const double e_i   = (t < 0.05) ? 0 : E_0 * std::sin(5 * (2 * M_PI) * (t - 0.05));
        const double e_i   = (t < 0.05) ? 0 : E_0;
        const double Tau_L = (t < 0.2) ? 0 : Tau_L0;
        dxdt[0]            = -a * x[0] - b * x[1] + e_i / L_A;
        dxdt[1]            = c * x[0] - d * x[1] - e * sign(x[1]) - Tau_L / Jeq;
        dxdt[2]            = f * x[0] * x[0] - g * x[2] + g * T_Amb;
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            current.emplace_back(x[0]);
            speed.emplace_back(x[1]);
            temperature.emplace_back(x[2]);
        }
    }
    std::vector<Variable> time, current, speed, temperature;
};

} // namespace dcmotor_v2

int main(int, char **)
{
    using namespace dcmotor_v2;

    // Instantiate the model.
    Model model;

    // Initial and runtime states.
    State x0{ .0, .0, 18.0 }, x;

    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    const auto Iterations = 10;
    using Euler           = chainsaw::stepper_euler<State, Time>;
    using Rk4             = chainsaw::stepper_rk4<State, Time>;
    using AdaptiveEuler   = chainsaw::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = chainsaw::stepper_adaptive<Rk4, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveEuler adaptive_euler;
    adaptive_euler.set_tollerance(1e-03);
    adaptive_euler.set_min_delta(1e-09);
    adaptive_euler.set_max_delta(1e-03);
    AdaptiveRk4 adaptive_rk4;
    adaptive_rk4.set_tollerance(1e-03);
    adaptive_rk4.set_min_delta(1e-09);
    adaptive_rk4.set_max_delta(1e-03);
    Euler euler;
    Rk4 rk4;

    // Setup the observers.
#ifdef ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs_adaptive_euler;
    Observer obs_adaptive_rk4;
    Observer obs_euler;
    Observer obs_rk4;

    // Instantiate the stopwatch.
    timelib::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << adaptive_euler.steps() << " steps, for a total of " << sw.partials()[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw.partials()[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << euler.steps() << " steps, for a total of " << sw.partials()[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << rk4.steps() << " steps, for a total of " << sw.partials()[3] << "\n";

#ifdef ENABLE_PLOT
    // Create a Gnuplot instance.
    gpcpp::Gnuplot gnuplot;
    
    // Set up the plot with grid, labels, and line widths
    gnuplot.set_title("Adaptive Methods vs Euler and RK4")
        .set_terminal(gpcpp::terminal_type_t::wxt)
        .set_xlabel("Time (s)")
        .set_ylabel("Values")
        .set_grid()
        .set_legend();

    // Plot scatter for Adaptive Euler - Current
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)       // Points style
        .set_point_style(gpcpp::point_style_t::filled_circle) // Marker style: filled circle ("o")
        .set_point_size(1)                             // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.current, "AdaptiveEuler.Current");

    // Plot scatter for Adaptive Euler - Speed
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)        // Points style
        .set_point_style(gpcpp::point_style_t::filled_diamond) // Marker style: filled circle ("o")
        .set_point_size(1)                              // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.speed, "AdaptiveEuler.Speed");

    // Plot scatter for Adaptive Euler - Temperature
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)       // Points style
        .set_point_style(gpcpp::point_style_t::filled_square) // Marker style: filled circle ("o")
        .set_point_size(1)                             // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.temperature, "AdaptiveEuler.Temperature");

    // Plot scatter for Adaptive RK4 - Current
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)     // Points style
        .set_point_style(gpcpp::point_style_t::open_circle) // Marker style: filled diamond ("d")
        .set_point_size(2)                           // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.current, "AdaptiveRk4.Current");

    // Plot scatter for Adaptive RK4 - Speed
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)      // Points style
        .set_point_style(gpcpp::point_style_t::open_diamond) // Marker style: filled diamond ("d")
        .set_point_size(2)                            // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.speed, "AdaptiveRk4.Speed");

    // Plot scatter for Adaptive RK4 - Temperature
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)     // Points style
        .set_point_style(gpcpp::point_style_t::open_square) // Marker style: filled diamond ("d")
        .set_point_size(2)                           // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.temperature, "AdaptiveRk4.Temperature");

    // Plot Euler method - Current
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_euler.time, obs_euler.current, "Euler.Current");

    // Plot Euler method - Speed
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_euler.time, obs_euler.speed, "Euler.Speed");

    // Plot Euler method - Temperature
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_euler.time, obs_euler.temperature, "Euler.Temperature");

    // Plot RK4 method - Current
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_rk4.time, obs_rk4.current, "Rk4.Current");

    // Plot RK4 method - Speed
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_rk4.time, obs_rk4.speed, "Rk4.Speed");

    // Plot RK4 method - Temperature
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_rk4.time, obs_rk4.temperature, "Rk4.Temperature");

    // Enable legend and display
    gnuplot.show();
#endif
    return 0;
}