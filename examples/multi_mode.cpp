/// @file multi_mode.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <exception>
#include <iomanip>
#include <iostream>
#include <timelib/stopwatch.hpp>

#ifdef ENABLE_PLOT
#include <gpcpp/gnuplot.hpp>
#endif

#include "defines.hpp"

#include <numint/detail/observer.hpp>
#include <numint/solver.hpp>
#include <numint/stepper/stepper_adaptive.hpp>
#include <numint/stepper/stepper_improved_euler.hpp>

namespace multi_mode
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
        : V(9.6)
        , R(8.4)
        , L(0.0084)
        , J(0.01)
        , Kd(0.25)
        , Ke(0.1785)
        , Kt(141.6 * Ke)
        , Fd(0.064)
        , Fs(0.035)
        , Ts(1)
        , Gr(_Gr)
        , R_Th(2.2)
        , C_Th(9 / R_Th)
        , T_Amb(22)
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
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        (void)t;
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
struct ObserverSave : public numint::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            current.emplace_back(x[0]);
            speed.emplace_back(x[1]);
            depth.emplace_back(x[2]);
            temperature.emplace_back(x[3]);
        }
    }

    std::vector<Variable> time, current, speed, depth, temperature;
};

} // namespace multi_mode

int main(int, char **)
{
    using namespace multi_mode;

    // Instantiate the model.
    Model model;
    // Initial and runtime states.
    State x0{.0, .0, .0, 18.0}, x;
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    // Setup the solvers.
    const auto Error      = numint::ErrorFormula::Mixed;
    const auto Iterations = 2;
    using Stepper         = numint::stepper_adaptive<numint::stepper_improved_euler<State, Time>, Iterations, Error>;
    // Instantiate the solvers.
    Stepper stepper;
    stepper.set_tollerance(1e-02);
    stepper.set_min_delta(1e-09);
    stepper.set_max_delta(1e-03);
    // Instantiate the observers.
#ifdef ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = numint::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer observer;
    // Instantiate the stopwatch.
    timelib::Stopwatch sw;
    x = x0;
    sw.start();
    numint::integrate_adaptive(stepper, observer, model, x, time_start, time_end / 2, time_delta);
    model.Gr = 10;
    numint::integrate_adaptive(stepper, observer, model, x, time_end / 2, time_end, time_delta);
    sw.round();
    std::cout << "Integration took " << std::setw(12) << stepper.steps() << " steps, for a total of " << sw.last_round()
              << "\n";

#ifdef ENABLE_PLOT
    // Create a Gnuplot instance.
    gpcpp::Gnuplot gnuplot;

    // Set up the plot with grid, labels, and line widths
    gnuplot.set_title("Observer Data")
        .set_terminal(gpcpp::terminal_type_t::wxt)
        .set_xlabel("Time (s)")
        .set_ylabel("Values")
        .set_grid()
        .set_legend();

    // Plot for Current
    gnuplot
        .set_line_width(2)                        // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .plot_xy(observer.time, observer.current, "Current");

    // Plot for Speed
    gnuplot
        .set_line_width(2)                        // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .plot_xy(observer.time, observer.speed, "Speed");

    // Plot for Temperature
    gnuplot
        .set_line_width(2)                        // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .plot_xy(observer.time, observer.temperature, "Temperature");

    gnuplot.show();

#endif
    return 0;
}