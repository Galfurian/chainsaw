/// @file model.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-04-13

#pragma once

#include "solver/observer.hpp"
#include "defines.hpp"

#include <exception>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>

namespace dcmotor_l
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Depth
/// x[3] : Temperature
using State = std::array<double, 4>;

/// @brief This one just containts the parameters.
struct Parameters {
    /// Supplied voltage[V].
    double V;
    /// Winding resistance in Ohms.
    double R;
    /// Winding inductance in Henrys[H].
    double L;
    /// Angular momentum[kg.m ^ 2].
    double J;
    /// Coulomb friction[N.m].
    double Kd;
    /// Back - EMF contanst[V * s / rad].
    double Ke;
    /// Torque constant[N * m / A].
    double Kt;
    /// Dynamic hole friction[Nm / mm]
    double Fd;
    /// Static hole  friction[Nm]
    double Fs;
    /// Thread slope, i.e., y - axis depth per revolution[mm / rev].
    double Ts;
    /// Gear ratio.
    double Gr;
    /// Thermal resistance of the motor [C / Watt].
    double R_Th;
    /// Thermal capacity of the coil [Joule / C].
    double C_Th;
    /// Ambient temperature.
    double T_Amb;

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

        ret.R_Th  = 9.9;
        ret.C_Th  = 133;
        ret.T_Amb = 22;
        return ret;
    }

    /// @brief Generates the parameters based on the given gear ratio.
    constexpr static Parameters params_n(double Gr) noexcept
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
    constexpr void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        /// x[0] : Current
        /// x[1] : Angular Speed
        /// x[2] : Depth
        /// x[3] : Temperature
        dxdt[0] = -(param.Kd / param.J) * x[0] + (param.Kt / param.J) * x[1] - ((param.Fd * param.Gr) / param.J) * x[2] - (param.Gr / param.J) * param.Fs;
        dxdt[1] = -(param.Ke / param.L) * x[0] - (param.R / param.L) * x[1] + (param.V / param.L);
        dxdt[2] = ((param.Ts * param.Gr) / (2 * M_PI)) * x[0];
        dxdt[3] = (param.R / param.C_Th) * x[0] * x[0] + (param.T_Amb - x[3]) / (param.C_Th * param.R_Th);
    }
};

/// @brief The dc motor itself.
template <int DECIMATION>
struct ObserverSave : public DecimationObserver<DECIMATION> {
    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> speed;
    std::vector<double> depth;
    std::vector<double> temperature;

    ObserverSave()
        : DecimationObserver<DECIMATION>(),
          time(),
          current(),
          speed(),
          depth(),
          temperature()
    {
        // Nothing to do.
    }

    void operator()(const State &x, const Time &t) noexcept
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
struct ObserverSafety {
    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> speed;
    std::vector<double> depth;
    std::vector<double> temperature;

    ObserverSafety()
        : time(),
          current(),
          speed(),
          depth(),
          temperature()
    {
        // Nothing to do.
    }

    void operator()(const State &x, const Time &t) noexcept
    {
        time.emplace_back(t);
        current.emplace_back(x[0]);
        speed.emplace_back(x[1]);
        depth.emplace_back(x[2]);
        temperature.emplace_back(x[3]);
    }
};

/// @brief The dc motor itself.
struct ObserverPrint {
    void operator()(const State &x, const Time &t)
    {
        std::cout << std::fixed << std::setprecision(4) << t << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << "\n";
    }
};

/// @brief The dc motor itself.
struct ObserverNone {
    void operator()(const State &x, const Time &t)
    {
        (void)x;
        (void)t;
    }
};

} // namespace dcmotor_l