/// @file compare_fixed.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <iostream>
#include <iomanip>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_improved_euler.hpp>
#include <chainsaw/stepper/stepper_trapezoidal.hpp>
#include <chainsaw/stepper/stepper_simpsons.hpp>
#include <chainsaw/stepper/stepper_midpoint.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

namespace comparison
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
(void) t;
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

} // namespace comparison

template <class Stepper, class System, class Observer>
inline void run_test_fixed_step(
    const std::string &name,
    Stepper &stepper,
    Observer &observer,
    System &&system,
    const typename Stepper::state_type &initial_state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type delta_time)
{
    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;
    // Initialize the state.
    typename Stepper::state_type x = initial_state;
    // Start the stopwatch.
    sw.start();
    // Run the integration.
    chainsaw::integrate_fixed(stepper, observer, system, x, start_time, end_time, delta_time);
    // Stop the stopwatch.
    sw.round();
    // Output the info.
    std::cout << "    " << std::setw(16) << name;
    std::cout << " took " << std::setw(12) << stepper.steps() << " steps,";
    std::cout << " for a total of " << sw.last_round() << "\n";
}

int main(int, char **)
{
    using namespace comparison;
    // Instantiate the model.
    Model model;
    // Initial and runtime states.
    State x0{ 10., 4. };
    // Simulation parameters.
    const Time start_time = 0.0, end_time = 2.0, delta_time = 0.05;

    // Instantiate the solvers and observers.
    chainsaw::stepper_euler<State, Time> euler;
    chainsaw::stepper_improved_euler<State, Time> improved_euler;
    chainsaw::stepper_midpoint<State, Time> midpoint;
    chainsaw::stepper_trapezoidal<State, Time> trapezoidal;
    chainsaw::stepper_simpsons<State, Time> simpsons;
    chainsaw::stepper_rk4<State, Time> rk4;
    chainsaw::stepper_rk4<State, Time> reference;

    // Setup the observers.
#ifdef SC_ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs_euler;
    Observer obs_improved_euler;
    Observer obs_midpoint;
    Observer obs_trapezoidal;
    Observer obs_simpsons;
    Observer obs_rk4;
    Observer obs_reference;

    // Run the integration.

    std::cout << "\n";
    std::cout << "Running integration...\n";
    run_test_fixed_step("euler", euler, obs_euler, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("improved_euler", improved_euler, obs_improved_euler, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("midpoint", midpoint, obs_midpoint, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("trapezoidal", trapezoidal, obs_trapezoidal, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("simpsons", simpsons, obs_simpsons, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("rk4", rk4, obs_rk4, model, x0, start_time, end_time, delta_time);
    run_test_fixed_step("reference", reference, obs_reference, model, x0, start_time, end_time, 1e-03);

#ifdef SC_ENABLE_PLOT
    matplot::hold(matplot::on);
    matplot::plot(obs_euler.time, obs_euler.x0)->line_width(2).line_style(":").display_name("euler.x0");
    matplot::plot(obs_improved_euler.time, obs_improved_euler.x0)->line_width(2).line_style("*").display_name("improved_euler.x0");
    matplot::plot(obs_midpoint.time, obs_midpoint.x0)->line_width(2).line_style("+").display_name("midpoint.x0");
    matplot::plot(obs_trapezoidal.time, obs_trapezoidal.x0)->line_width(2).line_style("-.").display_name("trapezoidal.x0");
    matplot::plot(obs_simpsons.time, obs_simpsons.x0)->line_width(3).line_style("-.").display_name("simpsons.x0");
    matplot::plot(obs_rk4.time, obs_rk4.x0)->line_width(2).line_style("--").display_name("rk4.x0");
    matplot::plot(obs_reference.time, obs_reference.x0)->line_width(2).line_style("--").display_name("reference.x0");
    matplot::legend(matplot::on);
    matplot::show();
#endif
    return 0;
}
