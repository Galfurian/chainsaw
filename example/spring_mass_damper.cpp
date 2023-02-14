/// @file spring_mass_damper.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "solver/detail/observer.hpp"
#include "solver/stepper/stepper_adaptive.hpp"
#include "solver/stepper/stepper_euler.hpp"
#include "solver/stepper/stepper_rk4.hpp"
#include "solver/solver.hpp"
#include "defines.hpp"

namespace spring_mass_damper
{

/// @brief State of the system.
///     x[0] : Position.
///     x[1] : Velocity.
using State = std::array<Variable, 2>;

/// @brief Parameters of our model.
struct Parameter {
    /// @brief Mass [kg].
    Variable m = 5.0;
    /// @brief Spring stiffness [N/m].
    Variable k = 40.0;
    /// @brief Damping constant.
    Variable c = 5;
};

struct Model : public Parameter {
    Model(Parameter parameter = Parameter())
        : Parameter(parameter)
    {
        // Nothing to do.
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        dxdt[0] = x[1];
        dxdt[1] = -c / m * x[1] - k / m * x[0];
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, p, v;
    constexpr inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            p.emplace_back(x[0]);
            v.emplace_back(x[1]);
        }
    }
};

} // namespace spring_mass_damper

int main(int, char **)
{
    using namespace spring_mass_damper;
    // Instantiate the model.
    Model model;
    // Change model's parameters.
    model.m = 1;
    model.c = 1;
    model.k = 1;
    // Runtime state.
    State x_f;
    State x_a;
    // Initial states.
    const State x0{ 1.0, 0.0 };
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 10;
    const Time time_delta = 1e-03;
    // Setup the fixed solver.
    using FixedSolver = solver::stepper_rk4<State, Time>;
    // Setup the adaptive solver.
    const auto Iterations = 3;
    const auto Error      = solver::ErrorFormula::Mixed;
    using AdaptiveSolver  = solver::stepper_adaptive<solver::stepper_rk4<State, Time>, Iterations, Error>;
    // Instantiate the solvers.
    FixedSolver solver_f;
    AdaptiveSolver solver_a;
    solver_a.set_tollerance(1e-09);
    solver_a.set_min_delta(1e-12);
    solver_a.set_max_delta(1e-01);
    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    ObserverSave<0> obs_f;
    ObserverSave<0> obs_a;
#elif 1
    solver::detail::ObserverPrint<0> obs_f;
    solver::detail::ObserverPrint<0> obs_a;
#endif
    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;
    std::cout << std::fixed;
    std::cout << "Simulating...\n";
    // Set the initial state.
    x_f = x0;
    x_a = x0;
    // Start the simulation.
    sw.start();
    // Run the solver.
    solver::integrate_fixed(solver_f, obs_f, model, x_f, time_start, time_end, time_delta);
    // Get the elapsed time.
    sw.round();
    // Run the solver.
    solver::integrate_adaptive(solver_a, obs_a, model, x_a, time_start, time_end, time_delta);
    // Get the elapsed time.
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Fixed solver computed    " << std::setw(12) << solver_f.steps() << " steps, for a total of " << sw[0] << "\n";
    std::cout << "    Adaptive solver computed " << std::setw(12) << solver_a.steps() << " steps, for a total of " << sw[1] << "\n";

#ifdef SC_ENABLE_PLOT
    auto figure = matplot::figure(true);
    matplot::grid(matplot::on);
    matplot::hold(matplot::on);
    matplot::plot(obs_f.time, obs_f.p)->line_width(2).display_name("Position F (m)");
    matplot::plot(obs_a.time, obs_a.p)->line_width(2).display_name("Position A (m)");
    matplot::plot(obs_f.time, obs_f.v, "--")->line_width(1).display_name("Speed F (m/s)");
    matplot::plot(obs_a.time, obs_a.v, "--")->line_width(1).display_name("Speed A (m/s)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}