#include <stopwatch/stopwatch.hpp>

#include "model.hpp"

#include "solver/solver.hpp"

#include <iostream>
#include <iomanip>

#include <matplot/matplot.h>

inline void run_fixed()
{
    using namespace dcmotor_l;
    auto param            = Parameters::default_params();
    auto dcmotor          = Model(param);
    auto x                = State{ .0, .0, .0, 22.0 };
    const Time time_start = 0.0;
    const Time time_end   = 100.0;
    const Time time_delta = 0.00001;
    const int samples     = compute_samples<int>(time_start, time_end, time_delta);

    solver::stepper_euler<State, Time> stepper;
    ObserverNone observer;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";
    std::cout << "Starting simulation.\n";

    stopwatch::StopWatch sw;
    auto steps = solver::integrate_const(stepper, observer, dcmotor, x, time_start, time_end, time_delta);
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n";
}

inline void run_adaptive()
{
    using namespace dcmotor_l;
    auto param             = Parameters::default_params();
    auto dcmotor           = Model(param);
    auto x                 = State{ .0, .0, .0, 22.0 };
    const Time time_start  = 0.0;
    const Time time_end    = 3.0;
    const Time time_delta  = 1e-16;
    const auto samples     = compute_samples<std::size_t>(time_start, time_end, time_delta);
    const auto downsamples = compute_samples<std::size_t>(time_start, time_end, time_delta);

    //solver::stepper_adaptive_euler<State, Time> stepper(0.00001);
    solver::stepper_adaptive_rk4<State, Time> stepper(1e-12);
    //ObserverNone observer;
    //ObserverPrint<0> observer;
    ObserverSave<0> observer;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";
    std::cout << "Starting simulation.\n";

    stopwatch::StopWatch sw;
    auto steps = solver::integrate_adaptive(stepper, observer, dcmotor, x, time_start, time_end, time_delta);
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n";

    matplot::hold(matplot::on);
    matplot::plot(observer.time, observer.current);
    matplot::plot(observer.time, observer.speed);
    matplot::plot(observer.time, observer.temperature);
    matplot::legend({ "Current", "Speed", "Temperature" });
    matplot::show();
}

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
    inline constexpr void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        dxdt[0] = 1.5 * x[0] - 1 * x[0] * x[1];
        dxdt[1] = -3 * x[1] + 1 * x[0] * x[1];
    }
};

template <std::size_t DECIMATION>
struct ObserverSave : public DecimationObserver<DECIMATION> {
    std::vector<Variable> time, x0, x1;
    ObserverSave() = default;
    inline constexpr void operator()(const State &x, const Time &t) noexcept
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

int main(int argc, char **argv)
{
    lotka::Model model;
    lotka::State x{ 10., 4. };
    const Time time_start  = 0.0;
    const Time time_end    = 100.0;
    const Time time_delta  = 0.00001;
    const auto samples     = compute_samples<std::size_t>(time_start, time_end, time_delta);
    const auto downsamples = compute_samples<std::size_t>(time_start, time_end, time_delta, 0.001);
    unsigned steps         = 0;

    solver::stepper_adaptive_euler<lotka::State, Time> adaptive_euler(0.00001);
    solver::stepper_adaptive_rk4<lotka::State, Time> adaptive_rk4(1e-12);
    solver::stepper_euler<lotka::State, Time> euler;
    solver::stepper_rk4<lotka::State, Time> rk4;

    //ObserverNone observer;
    //lotka::ObserverPrint<0> observer;
    lotka::ObserverSave<10> observer_adaptive_euler;
    lotka::ObserverSave<10> observer_adaptive_rk4;
    lotka::ObserverSave<downsamples> observer_euler;
    lotka::ObserverSave<downsamples> observer_rk4;

    stopwatch::StopWatch sw;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";

    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_euler, observer_adaptive_euler, model, x, time_start, time_end, 1e-12);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_rk4, observer_adaptive_rk4, model, x, time_start, time_end, 1e-12);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_const(euler, observer_euler, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = lotka::State{ 10., 4. };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_const(rk4, observer_rk4, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    auto colors      = matplot::palette::accent(4);
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

    //lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x1);
    //lh->line_width(3);
    //lh->color(matplot::to_array(colors[color_index++]));
    //lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x1);
    //lh->line_width(3);
    //lh->color(matplot::to_array(colors[color_index++]));
    //lh = matplot::plot(observer_euler.time, observer_euler.x1);
    //lh->line_width(3);
    //lh->color(matplot::to_array(colors[color_index++]));
    //lh = matplot::plot(observer_rk4.time, observer_rk4.x1);
    //lh->line_width(3);
    //lh->color(matplot::to_array(colors[color_index++]));

    matplot::legend(
        { 
            "adaptive_euler.x0",
            "adaptive_rk4.x0",
            "euler.x0",
            "rk4.x0",
            //"adaptive_euler.x1",
            //"adaptive_rk4.x1",
            //"euler.x1",
            //"rk4.x1" 
        }
    );
    matplot::show();

    //run_fixed();
    //run_adaptive();

    return 0;
}