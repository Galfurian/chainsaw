/// @file bounching_ball.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <timelib/stopwatch.hpp>
#include <cmath>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

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
    explicit Model(Parameter parameter)
        : Parameter(parameter)
    {
        // Nothing to do.
    }

    /// @brief Bouncing ball behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        (void)t;
        // x[0] : Ball velocity.
        // x[1] : Ball displacement.
        // Compute how much the ball has penetrated the ground.
        Variable penetration = x[1] - r;
        // Check for collision.
        if (penetration > 0) {
            dxdt[0] = -g;
        } else {
            dxdt[0] = -g - ((k * penetration) / m) - ((c * x[0]) / m);
        }
        dxdt[1] = x[0];
    }
};

template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            velocity.emplace_back(x[0]);
            displacement.emplace_back(x[1]);
        }
    }
    std::vector<Variable> time, velocity, displacement;
};

} // namespace bounching_ball

int main(int, char **)
{
    using namespace bounching_ball;
    //
    Parameter parameter;
    // Instantiate the model.
    Model model(parameter);
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
    const Time time_start = 0.0, time_end = 0.5, time_delta = 1e-02;
    // Setup the fixed solver.
    using FixedSolver = chainsaw::stepper_rk4<State, Time>;
    // Setup the adaptive solver.
    const auto Iterations = 16;
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    using AdaptiveSolver  = chainsaw::stepper_adaptive<FixedSolver, Iterations, Error>;

    // Instantiate the solvers.
    FixedSolver solver_f;
    AdaptiveSolver solver_a;
    solver_a.set_tollerance(1e-09);
    solver_a.set_min_delta(1e-12);
    solver_a.set_max_delta(1e-02);

    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    using Observer = ObserverSave<0>;
#elif 1
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs_f;
    Observer obs_a;

    // Instantiate the stopwatch.
    timelib::Stopwatch sw;
    std::cout << std::fixed;
    std::cout << "Simulating...\n";
    // Set the initial state.
    x_f = x0;
    x_a = x0;
    // Start the simulation.
    sw.start();
    // Run the solver.
    chainsaw::integrate_fixed(solver_f, obs_f, model, x_f, time_start, time_end, time_delta);
    // Get the elapsed time.
    sw.round();
    // Run the solver.
    chainsaw::integrate_adaptive(solver_a, obs_a, model, x_a, time_start, time_end, time_delta);
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
    matplot::line(0, 0, time_end, 0)->line_width(2).display_name("Ground");
    matplot::line(0, model.r, time_end, model.r)->line_width(1).line_style("--").display_name("Ball radius (m)");
    matplot::plot(obs_f.time, obs_f.displacement)->line_width(2).display_name("Position F (m)");
    matplot::plot(obs_a.time, obs_a.displacement)->line_width(2).display_name("Position A (m)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}