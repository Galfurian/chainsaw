/// @file dcmotor.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <iostream>
#include <iomanip>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "solver/detail/observer.hpp"
#include "solver/stepper/stepper_adaptive.hpp"
#include "solver/stepper/stepper_euler.hpp"
#include "solver/stepper/stepper_rk4.hpp"
#include "solver/solver.hpp"

#include "defines.hpp"

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
        ret.R  = 8.4;
        ret.L  = 0.0084;
        ret.J  = 0.01;
        ret.Kd = 0.25;
        ret.Ke = 0.1785;
        ret.Kt = 141.6 * ret.Ke;
        ret.Fd = 0.064;
        ret.Fs = 0.035;
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
    constexpr inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        /// x[0] : Current
        /// x[1] : Angular Speed
        /// x[2] : Depth
        /// x[3] : Temperature
        dxdt[0] = -(param.Kd / param.J) * x[0] + (param.Kt / param.J) * x[1] - ((param.Fd * param.Gr) / param.J) * x[2] - (param.Gr / param.J) * param.Fs;
        dxdt[1] = -(param.Ke / param.L) * x[0] - (param.R / param.L) * x[1] + (param.V / param.L);
        dxdt[2] = ((param.Ts * param.Gr) / (2 * M_PI)) * x[0];
        dxdt[3] = +(param.R / param.C_Th) * x[1] * x[1] + (param.T_Amb - x[3]) / (param.C_Th * param.R_Th);
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, current, speed, depth, temperature;

    ObserverSave() = default;

    constexpr inline void operator()(const State &x, const Time &t) noexcept
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
struct ObserverPrint : public solver::detail::DecimationObserver<DECIMATION> {
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
    using state_type_t = dcmotor::State;

    auto param = dcmotor::Parameters::default_params();
    dcmotor::Model model(param);
    state_type_t x0{ .0, .0, .0, 22.0 }, x;
    const Time time_start = 0.0;
#ifdef SC_ENABLE_PLOT
    const Time time_end = 1.0;
#else
    const Time time_end = 100.0;
#endif
    const Time time_delta  = 0.0001;
    const auto samples     = compute_samples<std::size_t>(time_start, time_end, time_delta);

    solver::stepper_adaptive<state_type_t, Time, solver::stepper_euler<state_type_t, Time>, 2> adaptive_euler(time_delta);
    solver::stepper_adaptive<state_type_t, Time, solver::stepper_rk4<state_type_t, Time>, 2> adaptive_rk4(time_delta);
    solver::stepper_euler<state_type_t, Time> euler;
    solver::stepper_rk4<state_type_t, Time> rk4;
    
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
    auto colors = matplot::palette::accent(16);
    auto color  = colors.begin();
    matplot::hold(matplot::on);
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.current)->line_width(3).color(matplot::to_array(*color++)).marker(matplot::line_spec::marker_style::asterisk);
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.current)->line_width(3).color(matplot::to_array(*color++)).marker(matplot::line_spec::marker_style::circle);
    matplot::plot(obs_euler.time, obs_euler.current)->line_width(3).color(matplot::to_array(*color++));
    matplot::plot(obs_rk4.time, obs_rk4.current)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.speed)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.speed)->line_width(3).color(matplot::to_array(*color++));
    matplot::plot(obs_euler.time, obs_euler.speed)->line_width(3).color(matplot::to_array(*color++));
    matplot::plot(obs_rk4.time, obs_rk4.speed)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.temperature)->line_width(3).color(matplot::to_array(*color++));
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.temperature)->line_width(3).color(matplot::to_array(*color++));
    matplot::plot(obs_euler.time, obs_euler.temperature)->line_width(3).color(matplot::to_array(*color++));
    matplot::plot(obs_rk4.time, obs_rk4.temperature)->line_width(3).color(matplot::to_array(*color++));
    // matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.depth)->line_width(3).color(matplot::to_array(*color++));
    // matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.depth)->line_width(3).color(matplot::to_array(*color++));
    // matplot::plot(obs_euler.time, obs_euler.depth)->line_width(3).color(matplot::to_array(*color++));
    // matplot::plot(obs_rk4.time, obs_rk4.depth)->line_width(3).color(matplot::to_array(*color++));

    matplot::legend(
        {
            "Adaptive Euler.current",
            "Adaptive RK4.current",
            "Euler.current",
            "RK4.current",
            "Adaptive Euler.speed",
            "Adaptive RK4.speed",
            "Euler.speed",
            "RK4.speed",
            "Adaptive Euler.temperature",
            "Adaptive RK4.temperature",
            "Euler.temperature",
            "RK4.temperature"
            // "Adaptive Euler.depth",
            // "Adaptive RK4.depth",
            // "Euler.depth",
            // "RK4.depth"
        });
    matplot::show();
#endif
}

int main(int, char **)
{
    compare_steppers();
    return 0;
}