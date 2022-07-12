#include <iomanip>
#include <iostream>

#include "solver/solver.hpp"
#include "model.hpp"

#include <stopwatch/stopwatch.hpp>
#include <matplot/matplot.h>

void run_dcmotor_plot()
{
    using namespace dcmotor_l;

    auto param            = Parameters::default_params();
    auto dcmotor          = Model(param);
    auto x                = State{ .0, .0, .0, 22.0 };
    const Time time_start = 0.0;
    const Time time_end   = 3.0;
    const Time time_delta = 0.0001;
    const int samples     = static_cast<int>(((time_end - time_start) / time_delta));
    const int downsamples = static_cast<int>(((time_end - time_start) / time_delta) * 0.01);

    //solver::stepper_euler<State, Time> stepper;
    solver::stepper_rk4<State, Time> stepper;

    //ObserverPrint observer;
    //ObserverNone observer;
    ObserverSave<downsamples> observer;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total time points " << samples << ".\n";
    std::cout << "Sampling with rate " << downsamples << ".\n";
    std::cout << "Starting simulation.\n";

    stopwatch::StopWatch sw;
    solver::integrate_const(stepper, observer, dcmotor, x, time_start, time_end, time_delta);
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";

    auto colors = matplot::palette::accent(4);

    matplot::line_handle lh;
    //matplot::hold(matplot::on);
    lh = matplot::plot(observer.time, observer.current);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[0]));
    matplot::ylabel("Current (A)");
    matplot::xlabel("Time (s)");
    matplot::save("current.png");

    lh = matplot::plot(observer.time, observer.speed);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[1]));
    matplot::ylabel("Speed (rad/s)");
    matplot::xlabel("Time (s)");
    matplot::save("speed.png");

    lh = matplot::plot(observer.time, observer.depth);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[2]));
    matplot::ylabel("Depth (mm)");
    matplot::xlabel("Time (s)");
    matplot::save("depth.png");

    solver::integrate_const(stepper, observer, dcmotor, x, time_end, 100, time_delta);

    lh = matplot::plot(observer.time, observer.temperature);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[3]));
    matplot::ylabel("Temperature (C)");
    matplot::xlabel("Time (s)");
    matplot::save("temperature.png");
}

void run_dcmotor()
{
    using namespace dcmotor_l;

    auto param            = Parameters::default_params();
    auto dcmotor          = Model(param);
    auto x                = State{ .0, .0, .0, 22.0 };
    const Time time_start = 0.0;
    const Time time_end   = 3.0;
    const Time time_delta = 0.0001;
    const int samples     = static_cast<int>(((time_end - time_start) / time_delta));
    const int downsamples = static_cast<int>(((time_end - time_start) / time_delta) * 0.01);

    //solver::stepper_euler<State, Time> stepper;
    solver::stepper_rk4<State, Time> stepper;

    //ObserverPrint<0> observer;
    ObserverNone observer;
    //ObserverSave<downsamples> observer;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Total time points " << samples << ".\n";
    std::cout << "Sampling with rate " << downsamples << ".\n";
    std::cout << "Starting simulation.\n";

    stopwatch::StopWatch sw;
    auto steps = solver::integrate_const(stepper, observer, dcmotor, x, time_start, time_end, time_delta);
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n";
}

void run_dcmotor_adaptive()
{
    using namespace dcmotor_l;

    auto param            = Parameters::default_params();
    auto dcmotor          = Model(param);
    auto x                = State{ .0, .0, .0, 22.0 };
    const Time time_start = 0.0;
    const Time time_end   = 3.0;
    const Time time_delta = 0.0001;
    const int samples     = static_cast<int>(((time_end - time_start) / time_delta));
    const int downsamples = static_cast<int>(((time_end - time_start) / time_delta) * 0.01);

    //solver::stepper_euler<State, Time> stepper;
    solver::stepper_adaptive<State, Time> stepper(0.00001, 0.01);

    //ObserverPrint<0> observer;
    //ObserverNone observer;
    ObserverSave<downsamples> observer;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << ".\n";
    std::cout << "Sampling with rate " << downsamples << ".\n";
    std::cout << "Starting simulation.\n";

    stopwatch::StopWatch sw;
    auto steps = solver::integrate_adaptive(stepper, observer, dcmotor, x, time_start, time_end, time_delta);
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n";
    
    auto colors = matplot::palette::accent(4);
    matplot::line_handle lh;
    //matplot::hold(matplot::on);
    lh = matplot::plot(observer.time, observer.current);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[0]));
    matplot::ylabel("Current (A)");
    matplot::xlabel("Time (s)");
    matplot::save("current.png");

    lh = matplot::plot(observer.time, observer.speed);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[1]));
    matplot::ylabel("Speed (rad/s)");
    matplot::xlabel("Time (s)");
    matplot::save("speed.png");

    lh = matplot::plot(observer.time, observer.depth);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[2]));
    matplot::ylabel("Depth (mm)");
    matplot::xlabel("Time (s)");
    matplot::save("depth.png");
}

int main(int argc, char **argv)
{
    std::cout << "\n";
    //run_dcmotor();
    run_dcmotor_plot();
    std::cout << "\n";
    run_dcmotor_adaptive();
    std::cout << "\n";
    return 0;
}