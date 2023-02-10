/// @file bounching_ball.cpp
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

namespace bounching_ball
{

/// @brief State of the system.
///     x[0] : Ball velocity.
///     x[1] : Ball displacement.
using State = std::array<Variable, 2>;

/// @brief Parameters of our model.
/// @details
/// The restitution coefficient is usually a positive, real number between 0 and 1:
///     e = 0     : Perfectly inelastic collision;
///     0 < e < 1 : Realistic inelastic collision, some kinetic energy is dissipated;
///     e = 1     : Full elastic collision, no kinetic energy is dissipated;
///     e > 1     : Energy is released uppon collision, e.g., nitrocellulose billiard balls;
struct Parameter {
    /// @brief Gravitational force [N].
    Variable g = 9.81;
    /// @brief Ball radius [m].
    Variable r = 0.01;
    /// @brief Ball mass [kg].
    Variable m = 0.2;
    /// @brief Spring constant.
    Variable k = 5000;
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
    constexpr inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        // Compute how much the ball has penetrated the ground.
        Variable penetration = x[1] + r;
        // Check for collision.
        if (penetration > 0) {
            dxdt[0] = -g;
        } else {
            dxdt[0] = -g - ((k * penetration) / m) - ((c * x[0]) / m);
        }
        dxdt[1] = x[0];
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, v, d;
    constexpr inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            v.emplace_back(x[0]);
            d.emplace_back(x[1]);
        }
    }
};

} // namespace bounching_ball

int main(int, char **)
{
    using namespace bounching_ball;
    // Instantiate the model.
    Model model;
    // Change model's parameters.
    model.g = 9.81;
    model.k = 5000;
    model.c = 5;
    model.m = 0.2;
    model.r = 0.01;
    // Runtime state.
    State x_f;
    State x_a;
    // Initial states.
    const State x0{ 0.0, 1.0 };
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 0.75;
    const Time time_delta = 1e-03;
    // Setup the fixed solver.
    using FixedSolver = solver::stepper_rk4<State, Time>;
    // Setup the adaptive solver.
    const auto Iterations = 10;
    const auto Error      = solver::ErrorFormula::Mixed;
    using AdaptiveSolver  = solver::stepper_adaptive<solver::stepper_rk4<State, Time>, Iterations, Error>;
    // Instantiate the solvers.
    FixedSolver solver_f;
    AdaptiveSolver solver_a;
    solver_a.set_tollerance(0.5);
    solver_a.set_min_delta(1e-12);
    solver_a.set_max_delta(1e-03);
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
    matplot::plot(obs_f.time, obs_f.d)->line_width(1).display_name("Position F (m)");
    matplot::plot(obs_a.time, obs_a.d)->line_width(1).display_name("Position A (m)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}