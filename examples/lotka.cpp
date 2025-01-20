/// @file lotka.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <timelib/stopwatch.hpp>
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

namespace lotka
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Depth
/// x[3] : Temperature
using State = std::array<Variable, 2>;

class Model {
public:
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        (void)t;
        dxdt[0] = 1.5 * x[0] - 1 * x[0] * x[1];
        dxdt[1] = -3 * x[1] + 1 * x[0] * x[1];
    }
};

template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            x0.emplace_back(x[0]);
            x1.emplace_back(x[1]);
        }
    }

    std::vector<Variable> time, x0, x1;
};

} // namespace lotka

int main(int, char **)
{
    using namespace lotka;

    // Instantiate the model.
    Model model;
    // Initial and runtime states.
    State x0{ 10., 4. }, x;

    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    const auto Iterations = 2;
    using Euler           = chainsaw::stepper_euler<State, Time>;
    using Rk4             = chainsaw::stepper_rk4<State, Time>;
    using AdaptiveEuler   = chainsaw::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = chainsaw::stepper_adaptive<Rk4, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveEuler adaptive_euler;
    adaptive_euler.set_tollerance(1e-03);
    adaptive_euler.set_min_delta(1e-09);
    adaptive_euler.set_max_delta(1e-02);
    AdaptiveRk4 adaptive_rk4;
    adaptive_rk4.set_tollerance(1e-03);
    adaptive_rk4.set_min_delta(1e-09);
    adaptive_rk4.set_max_delta(1e-02);
    Euler euler;
    Rk4 rk4;

    // Instantiate the observers.
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

    // Scatter plot for Adaptive Euler - x0
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)       // Points style
        .set_point_style(gpcpp::point_style_t::filled_circle) // Marker style: filled circle ("o")
        .set_point_size(1)                             // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.x0, "AdaptiveEuler.x0");

    // Scatter plot for Adaptive Euler - x1
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)       // Points style
        .set_point_style(gpcpp::point_style_t::filled_circle) // Marker style: filled circle ("o")
        .set_point_size(1)                             // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.x1, "AdaptiveEuler.x1");

    // Scatter plot for Adaptive RK4 - x0
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)      // Points style
        .set_point_style(gpcpp::point_style_t::open_diamond) // Marker style: open diamond ("d")
        .set_point_size(2)                            // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.x0, "AdaptiveRK4.x0");

    // Scatter plot for Adaptive RK4 - x1
    gnuplot.set_plot_style(gpcpp::plot_style_t::points)      // Points style
        .set_point_style(gpcpp::point_style_t::open_diamond) // Marker style: open diamond ("d")
        .set_point_size(2)                            // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.x1, "AdaptiveRK4.x1");

    // Line plot for Euler - x0
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_euler.time, obs_euler.x0, "Euler.x0");

    // Line plot for Euler - x1
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_euler.time, obs_euler.x1, "Euler.x1");

    // Line plot for RK4 - x0
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_rk4.time, obs_rk4.x0, "RK4.x0");

    // Line plot for RK4 - x1
    gnuplot.set_line_width(2)                // Line width
        .set_plot_style(gpcpp::plot_style_t::lines) // Line style
        .plot_xy(obs_rk4.time, obs_rk4.x1, "RK4.x1");

    gnuplot.show();

#endif
    return 0;
}