#include <iomanip>
#include <iostream>

#include "solver/solver.hpp"
#include "model.hpp"

#include <stopwatch/stopwatch.hpp>
#include <matplot/matplot.h>

void run_dcmotor_l()
{
    using namespace dcmotor_l;

    auto param            = Parameters::default_params();
    auto dcmotor          = Model(param);
    auto x                = State{ .0, .0, .0, 22.0 };
    const Time time_start = 0.0;
    const Time time_end   = 3.0;
    const Time dt         = 0.0001;
    const int samples     = static_cast<int>(((time_end - time_start) / dt));
    const int downsamples = static_cast<int>(((time_end - time_start) / dt) * 0.01);

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
    sw.set_print_mode(stopwatch::human);
    sw.start();
    for (Time t = time_start; t < time_end; t += dt) {
        stepper.do_step(dcmotor, x, t, dt);
        observer(x, t);
    }
    sw.stop();

    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    {
        matplot::plot(observer.time, observer.current);
        matplot::hold(matplot::on);
    }
    {
        matplot::plot(observer.time, observer.speed);
        matplot::hold(matplot::on);
    }
    {
        matplot::plot(observer.time, observer.temperature);
    }
    matplot::legend({ "Current", "Speed", "Power", "Temperature" });

    matplot::show();
}

int main(int argc, char **argv)
{
    //run_dcmotor_l_temp();
    run_dcmotor_l();
    return 0;
}