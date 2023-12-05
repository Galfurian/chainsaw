/// @file dcmotor.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <exception>
#include <iostream>
#include <iomanip>

#ifdef SC_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

#include "defines.hpp"

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

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
public:
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

    Parameters(Variable _Gr = 20)
        : V(9.6),
          R(8.4),
          L(0.0084),
          J(0.01),
          Kd(0.25),
          Ke(0.1785),
          Kt(141.6 * Ke),
          Fd(0.064),
          Fs(0.035),
          Ts(1),
          Gr(_Gr),
          R_Th(2.2),
          C_Th(9 / R_Th),
          T_Amb(22)
    {
        // Nothing to do.
    }
};

/// @brief The dc motor itself.
struct Model : public Parameters {
    Model(Parameters _param = Parameters())
        : Parameters(std::move(_param))
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
        dxdt[0] = -(R / L) * x[0] - (Ke / L) * x[1] + (V / L);
        dxdt[1] = +(Kt / J) * x[0] - (Kd / J) * x[1] - ((Fd * Gr) / J) * x[2] - (Gr / J) * Fs;
        dxdt[2] = ((Ts * Gr) / (2 * M_PI)) * x[1];
        dxdt[3] = +(R / C_Th) * x[0] * x[0] + (T_Amb - x[3]) / (C_Th * R_Th);
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::DecimationObserver<DECIMATION> {
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

} // namespace dcmotor

int main(int, char **)
{
    using namespace dcmotor;

    // Instantiate the model.
    Model model;

    // Initial and runtime states.
    State x0{ .0, .0, .0, 18.0 }, x;

    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    const auto Iterations = 2;
    using Euler           = chainsaw::stepper_euler<State, Time>;
    using Rk4             = chainsaw::stepper_rk4<State, Time>;
    using AdaptiveEuler   = chainsaw::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = chainsaw::stepper_adaptive<Rk4, Iterations, Error>;

    // Instantiate the solvers.
    AdaptiveEuler adaptive_euler;
    adaptive_euler.set_tollerance(1e-03);
    adaptive_euler.set_min_delta(1e-09);
    adaptive_euler.set_max_delta(1e-03);
    AdaptiveRk4 adaptive_rk4;
    adaptive_rk4.set_tollerance(1e-02);
    adaptive_rk4.set_min_delta(1e-09);
    adaptive_rk4.set_max_delta(1e-03);
    Euler euler;
    Rk4 rk4;

    // Instantiate the observers.
#ifdef SC_ENABLE_PLOT
    dcmotor::ObserverSave obs_adaptive_euler;
    dcmotor::ObserverSave obs_adaptive_rk4;
    dcmotor::ObserverSave obs_euler;
    dcmotor::ObserverSave obs_rk4;
#elif 1
    chainsaw::detail::NoObserver obs_adaptive_euler;
    chainsaw::detail::NoObserver obs_adaptive_rk4;
    chainsaw::detail::NoObserver obs_euler;
    chainsaw::detail::NoObserver obs_rk4;
#else
    chainsaw::detail::ObserverPrint<0> obs_adaptive_euler;
    chainsaw::detail::ObserverPrint<0> obs_adaptive_rk4;
    chainsaw::detail::ObserverPrint<0> obs_euler;
    chainsaw::detail::ObserverPrint<0> obs_rk4;
#endif

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << adaptive_euler.steps() << " steps, for a total of " << sw.partials()[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw.partials()[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << euler.steps() << " steps, for a total of " << sw.partials()[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << rk4.steps() << " steps, for a total of " << sw.partials()[3] << "\n";

#ifdef SC_ENABLE_PLOT
    matplot::hold(matplot::on);
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.current)->marker_size(13).marker_style("o").display_name("AdaptiveEuler.Current");
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.speed)->marker_size(13).marker_style("o").display_name("AdaptiveEuler.Speed");
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.temperature)->marker_size(13).marker_style("o").display_name("AdaptiveEuler.Temperature");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.current)->marker_size(26).marker_style("d").display_name("AdaptiveRk4.Current");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.speed)->marker_size(26).marker_style("d").display_name("AdaptiveRk4.Speed");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.temperature)->marker_size(26).marker_style("d").display_name("AdaptiveRk4.Temperature");
    matplot::plot(obs_euler.time, obs_euler.current)->line_width(2).display_name("Euler.Current");
    matplot::plot(obs_euler.time, obs_euler.speed)->line_width(2).display_name("Euler.Speed");
    matplot::plot(obs_euler.time, obs_euler.temperature)->line_width(2).display_name("Euler.Temperature");
    matplot::plot(obs_rk4.time, obs_rk4.current)->line_width(2).display_name("Rk4.Current");
    matplot::plot(obs_rk4.time, obs_rk4.speed)->line_width(2).display_name("Rk4.Speed");
    matplot::plot(obs_rk4.time, obs_rk4.temperature)->line_width(2).display_name("Rk4.Temperature");
    matplot::legend(matplot::on);
    matplot::show();
#endif
    return 0;
}