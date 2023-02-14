/// @file pendulum.cpp
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

namespace pendulum
{

/// @brief State of the system.
///     x[0] : Angle.
///     x[1] : Velocity.
using State = std::array<Variable, 2>;

/// @brief Parameters of our model.
struct Parameter {
    /// @brief Mass of the rod [kg].
    const Variable m;
    /// @brief Lenght of the rod [m].
    const Variable l;
    /// @brief Rotational damping [N.m].
    const Variable b;
    /// @brief Gravitational force [N].
    const Variable g;
    /// @brief Rod's moment of inertia about its center of mass.
    const Variable I;
    Parameter(Variable _m = 3.0,
              Variable _l = 0.19,
              Variable _b = 0.1,
              Variable _g = 9.81)
        : m(_m),
          l(_l),
          b(_b),
          g(_g),
          I((4. / 3.) * _m * _l * _l)
    {
        // Nothing to do.
    }
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
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
    {
#if 1
        const Variable u = (t < 3) ? 5 : 0;
#else
        const Variable u = .0;
#endif
        // Equation system.
        dxdt[0] = x[1];
        dxdt[1] = (u - m * g * l * x[0] - b * x[1]) / (I + m * l * l);
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, angle, velocity;
    constexpr inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            angle.emplace_back(x[0]);
            velocity.emplace_back(x[1]);
        }
    }
};

} // namespace pendulum

int main(int, char **)
{
    using namespace pendulum;
    // Instantiate the model.
    Model model;
    // Change model's parameters.
    // Runtime state.
    State x_f;
    State x_a;
    // Initial states.
    const State x0{ 0.0, 0.0 };
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 40;
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
    matplot::plot(obs_f.time, obs_f.angle)->line_width(2).display_name("Angle F (rad)");
    matplot::plot(obs_a.time, obs_a.angle)->line_width(2).display_name("Angle A (rad)");
    matplot::plot(obs_f.time, obs_f.velocity, "--")->line_width(1).display_name("Angular Speed F (rad/s)");
    matplot::plot(obs_a.time, obs_a.velocity, "--")->line_width(1).display_name("Angular Speed A (rad/s)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}