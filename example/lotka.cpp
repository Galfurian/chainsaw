/// @file lotka.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
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
    inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        dxdt[0] = 1.5 * x[0] - 1 * x[0] * x[1];
        dxdt[1] = -3 * x[1] + 1 * x[0] * x[1];
    }
};

template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, x0, x1;
    ObserverSave() = default;
    inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            x0.emplace_back(x[0]);
            x1.emplace_back(x[1]);
        }
    }
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
    const auto Error      = solver::ErrorFormula::Mixed;
    const auto Iterations = 2;
    using Euler           = solver::stepper_euler<State, Time>;
    using Rk4             = solver::stepper_rk4<State, Time>;
    using AdaptiveEuler   = solver::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = solver::stepper_adaptive<Rk4, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveEuler adaptive_euler;
    adaptive_euler.set_tollerance(1e-03);
    adaptive_euler.set_min_delta(1e-09);
    adaptive_euler.set_max_delta(1e-01);
    AdaptiveRk4 adaptive_rk4;
    adaptive_rk4.set_tollerance(1e-03);
    adaptive_rk4.set_min_delta(1e-09);
    adaptive_rk4.set_max_delta(1e-01);
    Euler euler;
    Rk4 rk4;

    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    lotka::ObserverSave obs_adaptive_euler;
    lotka::ObserverSave obs_adaptive_rk4;
    lotka::ObserverSave obs_euler;
    lotka::ObserverSave obs_rk4;
#else
    solver::detail::NoObserver obs_adaptive_euler;
    solver::detail::NoObserver obs_adaptive_rk4;
    solver::detail::NoObserver obs_euler;
    solver::detail::NoObserver obs_rk4;
#endif

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";
    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    solver::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    solver::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    solver::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    solver::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << adaptive_euler.steps() << " steps, for a total of " << sw.partials()[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw.partials()[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << euler.steps() << " steps, for a total of " << sw.partials()[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << rk4.steps() << " steps, for a total of " << sw.partials()[3] << "\n";

#ifdef SC_ENABLE_PLOT
    matplot::hold(matplot::on);
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x0)->marker_size(13).marker_style("o").display_name("AdaptiveEuler.x0");
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x1)->marker_size(13).marker_style("o").display_name("AdaptiveEuler.x1");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x0)->marker_size(26).marker_style("d").display_name("AdaptiveRK4.x0");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x1)->marker_size(26).marker_style("d").display_name("AdaptiveRK4.x1");
    matplot::plot(obs_euler.time, obs_euler.x0)->line_width(2).display_name("Euler.x0");
    matplot::plot(obs_euler.time, obs_euler.x1)->line_width(2).display_name("Euler.x1");
    matplot::plot(obs_rk4.time, obs_rk4.x0)->line_width(2).display_name("RK4.x0");
    matplot::plot(obs_rk4.time, obs_rk4.x1)->line_width(2).display_name("RK4.x1");
    matplot::legend(matplot::on);
    matplot::show();
#endif
    return 0;
}