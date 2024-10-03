/// @file accellerometer.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <cmath>
#include <stopwatch/stopwatch.hpp>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

namespace accellerometer
{

/// @brief State of the system.
///     Position
///     Velocity
using State = std::array<Variable, 2>;

/// @brief Parameters of our model.
struct Parameter {
    /// @brief Gravitational force [N].
    Variable g = 9.81;

    /// @brief seismic mass [Kg].
    Variable M = 0.16e-06;
    /// @brief spring stiffness [N/m].
    Variable K = 2.6455;
    /// @brief damping coefficient.
    Variable D = 4e-06;

    /// @brief spring stiffness [].
    Variable A = 220e-12;
    /// @brief spring stiffness [].
    Variable D0 = 1.5e-06;

    // /// @brief Mass [kg].
    // Variable m = 5.0;
    // /// @brief Spring stiffness [N/m].
    // Variable k = 40.0;
    // /// @brief Damping constant.
    // Variable c = 5;
};

struct Model : public Parameter {
    explicit Model(Parameter parameter)
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
        (void) t;
        double F = 4e-06 * sin(2.0 * M_PI * 1e04 * t);
        dxdt[0]  = x[1];
        dxdt[1]  = F; // * (1 / M) - x[0] * (K / M) - x[1] * (D / M);
    }
};

template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            position.emplace_back(x[0]);
            velocity.emplace_back(x[1]);
        }
    }
    std::vector<Variable> time, position, velocity;
};

} // namespace accellerometer

int main(int, char **)
{
    using accellerometer::Model;
    using accellerometer::ObserverSave;
    using accellerometer::Parameter;
    using accellerometer::State;
    // Instantiate the parameters.
    Parameter parameter;
    // Instantiate the model.
    Model model(parameter);
    // Runtime state.
    State x;
    // Initial states.
    const State x0{ 1.0, 0.0 };
    // Simulation parameters.
    const Time time_start = 0.0, time_end = 200e-06, time_delta = 1e-07;
    // Setup the base solver.
    using BaseSolver = chainsaw::stepper_rk4<State, Time>;
    // Instantiate the solvers.
    BaseSolver solver;
    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    using Observer = ObserverSave<0>;
#elif 1
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs;

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;
    std::cout << std::fixed;
    std::cout << "Simulating...\n";

    // Set the initial state.
    x = x0;

    // Start the simulation.
    sw.start();
    // Run the solver.
    chainsaw::integrate_fixed(solver, obs, model, x, time_start, time_end, time_delta);
    // Get the elapsed time.
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Fixed solver computed    " << std::setw(12) << solver.steps() << " steps, for a total of " << sw[0] << "\n";

#ifdef SC_ENABLE_PLOT
    auto figure = matplot::figure(true);
    matplot::grid(matplot::on);
    matplot::hold(matplot::on);
    matplot::line(0, 0, time_end, 0)->line_width(2).display_name("Ground");
    matplot::plot(obs.time, obs.position)->line_width(2).display_name("Position (m)");
    matplot::plot(obs.time, obs.velocity)->line_width(2).display_name("Velocity (m/s)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}