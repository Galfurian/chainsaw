/// @file tandem_dc_motors.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include <stopwatch/stopwatch.hpp>
#include <exception>
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

namespace tandem_dc_motors
{

/// @brief State of the system.
///     x[0] : Current motor 1.
///     x[1] : Current motor 2.
///     x[2] : Angular speed motor 1.
///     x[3] : Angular speed motor 2.
///     x[4] : Angular speed load.
///     x[5] : Angle shaft 1.
///     x[6] : Angle shaft 2.
using State = std::array<Variable, 7>;

/// @brief The model.
struct Model {
    Variable v_a1;
    Variable v_a2;

    Variable R_a1;
    Variable L_a1;
    Variable K_e1;

    Variable R_a2;
    Variable L_a2;
    Variable K_e2;

    Variable Kt_m1;
    Variable J_m1;
    Variable Kd_m1;

    Variable Kt_m2;
    Variable J_m2;
    Variable Kd_m2;

    Variable Kd_s1;
    Variable Ke_s1;

    Variable Kd_s2;
    Variable Ke_s2;

    Variable J_l;
    Variable Kd_l;

    Model()
        : v_a1(50),
          v_a2(0),

          R_a1(8.4),
          L_a1(0.0084),
          K_e1(0.1785),

          R_a2(8.4),
          L_a2(0.0084),
          K_e2(0.1785),

          Kt_m1(141.6 * K_e1),
          J_m1(0.5),
          Kd_m1(0.05),

          Kt_m2(141.6 * K_e2),
          J_m2(0.5),
          Kd_m2(0.05),

          Kd_s1(0.05),
          Ke_s1(0.01),

          Kd_s2(0.05),
          Ke_s2(0.01),

          J_l(0.5),
          Kd_l(0.05)
    {
        // Nothing to do.
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    constexpr inline void operator()(const State &x, State &dxdt, Time) noexcept
    {
        /// x[0] : Current motor 1.
        /// x[1] : Current motor 2.
        /// x[2] : Angular speed motor 1.
        /// x[3] : Angular speed motor 2.
        /// x[4] : Angular speed load.
        /// x[5] : Angle shaft 1.
        /// x[6] : Angle shaft 2.
        dxdt[0] = -R_a1 / L_a1 * x[0] - K_e1 / L_a1 * x[2] + v_a1;
        dxdt[1] = -R_a2 / L_a2 * x[1] - K_e2 / L_a2 * x[3] + v_a2;

        dxdt[2] = Kt_m1 / J_m1 * x[0] - (Kd_m1 + Kd_s1) / J_m1 * x[2] + Kd_s2 / J_m1 * x[3] - Ke_s1 / J_m1 * x[5];
        dxdt[3] = Kt_m2 / J_m2 * x[1] + Kd_s1 / J_m2 * x[2] - (Kd_m1 + Kd_m2 + Kd_s1) / J_m2 * x[3] + Kd_s2 / J_m2 * x[4] + Ke_s1 / J_m2 * x[5] - Ke_s2 / J_m2 * x[6];

        dxdt[4] = Kt_m1 / J_l * x[0] + Kt_m2 / J_l * x[1] + Kd_s2 / J_l * x[3] - (Kd_l + Kd_s2) / J_l * x[4] + Ke_s2 / J_l * x[6];
        
        dxdt[5] = x[2] - x[3];
        dxdt[6] = x[3] - x[4];
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, i_a1, i_a2, w_m1, w_m2, w_l, a_s1, a_s2;

    ObserverSave() = default;

    constexpr inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            i_a1.emplace_back(x[0]);
            i_a2.emplace_back(x[1]);
            w_m1.emplace_back(x[2]);
            w_m2.emplace_back(x[3]);
            w_l.emplace_back(x[4]);
            a_s1.emplace_back(x[5]);
            a_s2.emplace_back(x[6]);
        }
    }
};

/// @brief The dc motor itself.
template <std::size_t DECIMATION = 0>
struct ObserverPrint : public solver::detail::DecimationObserver<DECIMATION> {
    ObserverPrint() = default;
    inline void operator()(const State &x, const Time &t)
    {
        if (this->observe())
            std::cout << std::fixed << std::setprecision(4) << t << " " << x << "\n";
    }
};

} // namespace tandem_dc_motors

int main(int, char **)
{
    using namespace tandem_dc_motors;
    using state_type_t = State;

    Model model;
    state_type_t x0{ .0, .0, .0, .0, .0, .0, .0 }, x;
    const Time time_start = 0.0;
    const Time time_end   = 50.0;
    const Time time_delta = 0.000001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    using Rk4             = solver::stepper_rk4<state_type_t, Time>;
    const auto Error      = solver::ErrorFormula::Mixed;
    const auto Iterations = 2;

    solver::stepper_adaptive<state_type_t, Time, Rk4, Iterations, Error> adaptive_rk4(time_delta);

    std::size_t steps_adaptive_rk4;

#ifdef SC_ENABLE_PLOT
    ObserverSave obs_adaptive_rk4;
#elif 1
    ObserverPrint obs_adaptive_rk4;
#endif

    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    steps_adaptive_rk4 = solver::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << steps_adaptive_rk4 << " steps, for a total of " << sw.partials()[0] << "\n";

#ifdef SC_ENABLE_PLOT
    auto colors = matplot::palette::accent(7);
    auto color  = colors.begin();
    matplot::hold(matplot::on);

    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.i_a1)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.i_a2)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_m1)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_m2)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.w_l)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.a_s1)->color(matplot::to_array(*color++)).line_width(2.0);
    matplot::plot(obs_adaptive_rk4.time, obs_adaptive_rk4.a_s2)->color(matplot::to_array(*color++)).line_width(2.0);

    matplot::legend(
        { "Current M1",
          "Current M2",
          "Angular Speed M1",
          "Angular Speed M2",
          "Angular Speed Load",
          "Angle Shaft 1",
          "Angle Shaft 2" });
    matplot::show();
#endif
    return 0;
}