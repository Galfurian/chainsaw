/// @file lotka.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <matplot/matplot.h>

#include "solver/observer.hpp"
#include "solver/solver.hpp"
#include "defines.hpp"

#include <iostream>
#include <iomanip>

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
        dxdt[0] = 1.5 * x[0] - 1 * x[0] * x[1];
        dxdt[1] = -3 * x[1] + 1 * x[0] * x[1];
    }
};

template <std::size_t DECIMATION>
struct ObserverSave : public DecimationObserver<DECIMATION> {
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
struct ObserverPrint : public DecimationObserver<DECIMATION> {
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
    lotka::State x{ 10., 4. };
    const Time time_start  = 0.0;
    const Time time_end    = 10.0;
    const Time time_delta  = 0.00001;
    const auto samples     = compute_samples<std::size_t>(time_start, time_end, time_delta);
    const auto downsamples = compute_samples<std::size_t>(time_start, time_end, time_delta, 0.0001);
    unsigned steps         = 0;

    solver::stepper_adaptive_euler<lotka::State, Time> adaptive_euler(0.00001);
    solver::stepper_adaptive_rk4<lotka::State, Time, 2> adaptive_rk4(1e-12);
    solver::stepper_euler<lotka::State, Time> euler;
    solver::stepper_rk4<lotka::State, Time> rk4;

    //ObserverNone observer;
    //lotka::ObserverPrint<0> observer;
    lotka::ObserverSave<0> observer_adaptive_euler;
    lotka::ObserverSave<0> observer_adaptive_rk4;
    lotka::ObserverSave<downsamples> observer_euler;
    lotka::ObserverSave<downsamples> observer_rk4;

    stopwatch::Stopwatch sw;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";

    std::cout << "Starting simulation (Adaptive Euler).\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_euler, observer_adaptive_euler, model, x, time_start, time_end, 1e-16);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation (Adaptive RK4).\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_rk4, observer_adaptive_rk4, model, x, time_start, time_end, 1e-12);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation (Euler).\n";
    sw.start();
    steps = solver::integrate_const(euler, observer_euler, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation (RK4).\n";
    sw.start();
    steps = solver::integrate_const(rk4, observer_rk4, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    return;

    auto colors      = matplot::palette::accent(8);
    auto color_index = 0u;
    matplot::line_handle lh;
    matplot::hold(matplot::on);
    lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_euler.time, observer_euler.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_rk4.time, observer_rk4.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_euler.time, observer_euler.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_rk4.time, observer_rk4.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    matplot::legend(
        {
            "Adaptive Euler.x0",
            "Adaptive RK4.x0",
            "Euler.x0",
            "RK4.x0",
            "Adaptive Euler.x1",
            "Adaptive RK4.x1",
            "Euler.x1",
            "RK4.x1"
        });
    matplot::show();
}

int main(int argc, char **argv)
{
    compare_steppers();
    return 0;
}