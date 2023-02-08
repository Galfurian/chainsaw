/// @file dcmotor_nl_frict.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief
/// @version 0.1
/// @date 2022-04-13

#pragma once

#include "defines.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace dcmotor_nl_frict
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Depth
using State = std::array<double, 3>;

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
        lhs << "V   :" << rhs.V << "\n";
        lhs << "R   :" << rhs.R << "\n";
        lhs << "L   :" << rhs.L << "\n";
        lhs << "J   :" << rhs.J << "\n";
        lhs << "Kd  :" << rhs.Kd << "\n";
        lhs << "Ke  :" << rhs.Ke << "\n";
        lhs << "Kt  :" << rhs.Kt << "\n";
        lhs << "Fd  :" << rhs.Fd << "\n";
        lhs << "Fs  :" << rhs.Fs << "\n";
        lhs << "Ts  :" << rhs.Ts << "\n";
        lhs << "Gr  :" << rhs.Gr << "\n";
        return lhs;
    }

    constexpr inline double nonlinear_friction(const State &x) const
    {
        // nonlinear_friction Computes the friction in contact between moving bodies.
        // Inputs:
        //   Vr : is relative velocity between the two bodies.
        // Outputs:
        //   F  : is friction force.
        constexpr double Fc = 0.01;   // Coulomb friction.
        constexpr double Fb = 0.05;  // Breakaway friction.
        constexpr double Fv = 0.05; // Viscous friction.
        constexpr double Vs = 0.1; // Stribeck velocity.
        // Compute the friction.
        return (Fc + (Fb - Fc) * std::exp(-std::abs(x[0] / Vs))) * sign(x[0]) + Fv * x[0];
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
        //dxdt[0] = -(param.Kd / param.J) * x[0] + (param.Kt / param.J) * x[1] - ((param.Fd * param.Gr) / param.J) * x[2] - ((param.Fs * param.Gr) / param.J);
        dxdt[0] = -(param.Kd / param.J) * x[0] + (param.Kt / param.J) * x[1] - ((param.nonlinear_friction(x) * param.Gr) / param.J) * x[2];
        dxdt[1] = -(param.Ke / param.L) * x[0] - (param.R / param.L) * x[1] + (1 / param.L) * param.V;
        dxdt[2] = ((param.Ts * param.Gr) / (2 * M_PI)) * x[0];
    }
};

} // namespace dcmotor_nl_frict