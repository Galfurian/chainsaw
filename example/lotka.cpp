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
#include "solver/stepper/stepper_adaptive_euler.hpp"
#include "solver/stepper/stepper_adaptive_rk4.hpp"
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

template <std::size_t DECIMATION>
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

/// @brief The dc motor itself.
template <std::size_t DECIMATION>
struct ObserverPrint : public solver::detail::DecimationObserver<DECIMATION> {
    ObserverPrint() = default;
    inline void operator()(const State &x, const Time &t)
    {
        if (this->observe())
            std::cout << std::fixed << std::setprecision(4) << t << " " << x << "\n";
    }
};

} // namespace lotka

void compare_steppers()
{
    lotka::Model model;
    lotka::State x0{ 10., 4. }, x;
    const Time time_start = 0.0;
#ifdef SC_ENABLE_PLOT
    const Time time_end = 1.0;
#else
    const Time time_end = 100.0;
#endif
    const Time time_delta = 0.00001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    solver::stepper_adaptive_euler<lotka::State, Time> adaptive_euler(time_delta);
    solver::stepper_adaptive_rk4<lotka::State, Time> adaptive_rk4(time_delta);
    solver::stepper_euler<lotka::State, Time> euler;
    solver::stepper_rk4<lotka::State, Time> rk4;

    std::size_t steps_adaptive_euler;
    std::size_t steps_adaptive_rk4;
    std::size_t steps_euler;
    std::size_t steps_rk4;

#ifdef SC_ENABLE_PLOT
    dcmotor::ObserverSave obs_adaptive_euler;
    dcmotor::ObserverSave obs_adaptive_rk4;
    dcmotor::ObserverSave obs_euler;
    dcmotor::ObserverSave obs_rk4;
#else
    solver::detail::NoObserver obs_adaptive_euler;
    solver::detail::NoObserver obs_adaptive_rk4;
    solver::detail::NoObserver obs_euler;
    solver::detail::NoObserver obs_rk4;
#endif

    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    steps_adaptive_euler = solver::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    steps_adaptive_rk4 = solver::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    steps_euler = solver::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    steps_rk4 = solver::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << steps_adaptive_euler << " steps, for a total of " << sw.partials()[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << steps_adaptive_rk4 << " steps, for a total of " << sw.partials()[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << steps_euler << " steps, for a total of " << sw.partials()[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << steps_rk4 << " steps, for a total of " << sw.partials()[3] << "\n";

#ifdef SC_ENABLE_PLOT
    auto colors = matplot::palette::accent(8);
    auto color  = colors.begin();

    matplot::hold(matplot::on);

    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x0)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x0)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_euler.time, obs_euler.x0)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_rk4.time, obs_rk4.x0)->line_width(3).color(matplot::to_array(*color++));

    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x1)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x1)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_euler.time, obs_euler.x1)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_rk4.time, obs_rk4.x1)->line_width(3).color(matplot::to_array(*color++));

    matplot::legend(
        { "Adaptive Euler.x0",
          "Adaptive RK4.x0",
          "Euler.x0",
          "RK4.x0",
          "Adaptive Euler.x1",
          "Adaptive RK4.x1",
          "Euler.x1",
          "RK4.x1" });
    matplot::show();
#endif
}

int main(int, char **)
{
    compare_steppers();
    return 0;
}