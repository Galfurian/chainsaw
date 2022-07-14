/// @file dcmotor.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>

#include "solver/observer.hpp"
#include "solver/solver.hpp"
#include "defines.hpp"

#include <iostream>
#include <iomanip>

namespace dcmotor
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Depth
/// x[3] : Temperature
using State = std::array<Variable, 4>;

/// @brief This one just containts the parameters.
struct Parameters {
    /// Supplied voltage[V].
    Variable V;
    /// Winding resistance in Ohms.
    Variable R;
    /// Winding inductance in Henrys[H].
    Variable L;
    /// Angular momentum[kg.m ^ 2].
    Variable J;
    /// Coulomb friction[N.m].
    Variable Kd;
    /// Back - EMF contanst[V * s / rad].
    Variable Ke;
    /// Torque constant[N * m / A].
    Variable Kt;
    /// Dynamic hole friction[Nm / mm]
    Variable Fd;
    /// Static hole  friction[Nm]
    Variable Fs;
    /// Thread slope, i.e., y - axis depth per revolution[mm / rev].
    Variable Ts;
    /// Gear ratio.
    Variable Gr;
    /// Thermal resistance of the motor [C / Watt].
    Variable R_Th;
    /// Thermal capacity of the coil [Joule / C].
    Variable C_Th;
    /// Ambient temperature.
    Variable T_Amb;

    /// @brief Generates the default parameters.
    constexpr static Parameters default_params() noexcept
    {
        Parameters ret{};
        ret.V  = 9.6;
        ret.R  = 1;
        ret.L  = 25e-05;
        ret.J  = 0.1;
        ret.Kd = 0.25;
        ret.Ke = 1.00;
        ret.Kt = 1.00;
        ret.Fd = 0.01;
        ret.Fs = 0.05;
        ret.Ts = 1;
        ret.Gr = 20;

        ret.R_Th  = 2.2;
        ret.C_Th  = 9 / ret.R_Th;
        ret.T_Amb = 22;
        return ret;
    }

    /// @brief Generates the parameters based on the given gear ratio.
    constexpr static Parameters params_n(Variable Gr) noexcept
    {
        Parameters ret = Parameters::default_params();
        ret.Gr         = Gr;
        return ret;
    }

    friend std::ostream &operator<<(std::ostream &lhs, const Parameters &rhs)
    {
        lhs << "V       :" << rhs.V << "\n";
        lhs << "R       :" << rhs.R << "\n";
        lhs << "L       :" << rhs.L << "\n";
        lhs << "J       :" << rhs.J << "\n";
        lhs << "Kd      :" << rhs.Kd << "\n";
        lhs << "Ke      :" << rhs.Ke << "\n";
        lhs << "Kt      :" << rhs.Kt << "\n";
        lhs << "Fd      :" << rhs.Fd << "\n";
        lhs << "Fs      :" << rhs.Fs << "\n";
        lhs << "Ts      :" << rhs.Ts << "\n";
        lhs << "Gr      :" << rhs.Gr << "\n";
        lhs << "R_Th    :" << rhs.R_Th << "\n";
        lhs << "T_Amb   :" << rhs.T_Amb << "\n";
        return lhs;
    }
};

/// @brief The dc motor itself.
struct Model {
    /// The motor paramters.
    const Parameters param;

    Model(const Parameters &_param)
        : param(_param)
    {
        // Nothing to do.
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    inline constexpr void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        /// x[0] : Current
        /// x[1] : Angular Speed
        /// x[2] : Depth
        /// x[3] : Temperature
        dxdt[0] = -(param.Kd / param.J) * x[0] + (param.Kt / param.J) * x[1] - ((param.Fd * param.Gr) / param.J) * x[2] - (param.Gr / param.J) * param.Fs;
        dxdt[1] = -(param.Ke / param.L) * x[0] - (param.R / param.L) * x[1] + (param.V / param.L);
        dxdt[2] = ((param.Ts * param.Gr) / (2 * M_PI)) * x[0];
        dxdt[3] =
            +(param.R / param.C_Th) * x[1] * x[1] + (param.T_Amb - x[3]) / (param.C_Th * param.R_Th);
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION>
struct ObserverSave : public DecimationObserver<DECIMATION> {
    std::vector<Variable> time, current, speed, depth, temperature;

    ObserverSave() = default;

    inline constexpr void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            current.emplace_back(x[0]);
            speed.emplace_back(x[1]);
            depth.emplace_back(x[2]);
            temperature.emplace_back(x[3]);
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

} // namespace dcmotor

void compare_steppers()
{
    auto param = dcmotor::Parameters::default_params();
    dcmotor::Model model(param);
    dcmotor::State x;
    const Time time_start  = 0.0;
    const Time time_end    = 10.0;
    const Time time_delta  = 0.00001;
    const auto samples     = compute_samples<std::size_t>(time_start, time_end, time_delta);
    const auto downsamples = compute_samples<std::size_t>(time_start, time_end, time_delta, 0.001);
    unsigned steps         = 0;

    solver::stepper_adaptive_euler<dcmotor::State, Time> adaptive_euler(0.00001);
    solver::stepper_adaptive_rk4<dcmotor::State, Time> adaptive_rk4(1e-12);
    solver::stepper_euler<dcmotor::State, Time> euler;
    solver::stepper_rk4<dcmotor::State, Time> rk4;

    //ObserverNone observer;
    //dcmotor::ObserverPrint<0> observer;
    dcmotor::ObserverSave<10> observer_adaptive_euler;
    dcmotor::ObserverSave<10> observer_adaptive_rk4;
    dcmotor::ObserverSave<downsamples> observer_euler;
    dcmotor::ObserverSave<downsamples> observer_rk4;

    stopwatch::StopWatch sw;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";

    x = dcmotor::State{ .0, .0, .0, 22.0 };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_euler, observer_adaptive_euler, model, x, time_start, time_end, 1e-16);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = dcmotor::State{ .0, .0, .0, 22.0 };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_rk4, observer_adaptive_rk4, model, x, time_start, time_end, 1e-16);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = dcmotor::State{ .0, .0, .0, 22.0 };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_const(euler, observer_euler, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = dcmotor::State{ .0, .0, .0, 22.0 };
    std::cout << "Starting simulation.\n";
    sw.start();
    steps = solver::integrate_const(rk4, observer_rk4, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    // auto colors      = matplot::palette::accent(4);
    // auto color_index = 0u;
    // matplot::line_handle lh;
    // matplot::hold(matplot::on);
    // lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x0);
    // lh->line_width(3);
    // lh->color(matplot::to_array(colors[color_index++]));
    // lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x0);
    // lh->line_width(3);
    // lh->color(matplot::to_array(colors[color_index++]));
    // lh = matplot::plot(observer_euler.time, observer_euler.x0);
    // lh->line_width(3);
    // lh->color(matplot::to_array(colors[color_index++]));
    // lh = matplot::plot(observer_rk4.time, observer_rk4.x0);
    // lh->line_width(3);
    // lh->color(matplot::to_array(colors[color_index++]));
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
    // matplot::legend(
    //     {
    //         "adaptive_euler.x0",
    //         "adaptive_rk4.x0",
    //         "euler.x0",
    //         "rk4.x0",
    //         //"adaptive_euler.x1",
    //         //"adaptive_rk4.x1",
    //         //"euler.x1",
    //         //"rk4.x1"
    //     }
    // );
    // matplot::show();
    //run_fixed();
    //run_adaptive();
}

int main(int argc, char **argv)
{
    compare_steppers();
    return 0;
}