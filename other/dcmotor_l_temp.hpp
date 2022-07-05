/// @file dcmotor_l_temp.hpp
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

namespace dcmotor_l_temp
{

/// @brief State of the system.
/// x[0] : Current
/// x[1] : Angular Speed
/// x[2] : Depth
/// x[3] : Temperature
using State = std::array<double, 4>;

/// @brief This one just containts the parameters.
struct Parameters {
    double E_0;    // [V]
    double Tau_L0; // [N.m]
    double T_Amb;  // [deg]
    double B_2C;   // [N]

    // motor parameters , Nachtigal , Table 16.5 p. 663
    double J_1; // in*oz*s ^2/ rad
    double B_1; // in*oz*s/ rad

    // electrical / mechanical relations
    double K_E; // back emf coefficient , e_m = K_E * omega_m ( K_E= alpha * omega )
    double K_T; // torque coeffic ., in English units K_T is not = K_E ! ( K_T = alpha * )
    double R_A; // Ohms
    double L_A; // H

    // gear - train and load parameters
    double J_2; // in*oz*s^2/ rad // 10x motor J
    double B_2; // in*oz*s/ rad ( viscous )
    double N;   // motor / load gear ratio ; omega_1 = N omega_2

    // Thermal model parameters
    double R_TM; // Thermal resistance (C / Watt)
    double C_TM; // Thermal capacity (Watt - sec / C) (-> 9 sec time constant - fast !)

    // Support variables.
    double Jeq;
    double Beq;
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;

    /// @brief Generates the default parameters.
    constexpr static Parameters default_params() noexcept
    {
        Parameters param{};

        param.E_0    = 120.;
        param.Tau_L0 = 80.;
        param.T_Amb  = 18.;
        param.B_2C   = 300.;

        // motor parameters , Nachtigal , Table 16.5 p. 663
        param.J_1 = 0.0035;
        param.B_1 = 0.064;

        // electrical / mechanical relations
        param.K_E = 0.1785;
        param.K_T = 141.6 * param.K_E;
        param.R_A = 8.4;
        param.L_A = 0.0084;

        // gear - train and load parameters
        param.J_2 = 0.035;
        param.B_2 = 2.64;
        param.N   = 8.;

        // Thermal model parameters
        param.R_TM = 2.2;
        param.C_TM = 9. / param.R_TM;

        // Support variables.
        param.Jeq = param.J_2 + param.N * 2 * param.J_1;
        param.Beq = param.B_2 + param.N * param.N * param.B_1;
        param.a   = param.R_A / param.L_A;
        param.b   = param.K_E * param.N / param.L_A;
        param.c   = param.N * param.K_T / param.Jeq;
        param.d   = param.Beq / param.Jeq;
        param.e   = param.B_2C / param.Jeq;
        param.f   = param.R_A / param.C_TM;
        param.g   = 1. / (param.C_TM * param.R_TM);
        return param;
    }

    friend std::ostream &operator<<(std::ostream &lhs, const Parameters &rhs)
    {
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
        //const double e_i   = (t < 0.05) ? 0 : param.E_0 * std::sin(5 * (2 * M_PI) * (t - 0.05));
        const double e_i   = (t < 0.05) ? 0 : param.E_0;
        const double Tau_L = (t < 0.2) ? 0 : param.Tau_L0;
        /// x[0] : Current
        /// x[1] : Angular Speed
        /// x[2] : Temperature
        dxdt[0] = -param.a * x[0] - param.b * x[1] + e_i / param.L_A;
        dxdt[1] = param.c * x[0] - param.d * x[1] - param.e * sign(x[1]) - Tau_L / param.Jeq;
        dxdt[2] = param.f * x[0] * x[0] - param.g * x[2] + param.g * param.T_Amb;
    }
};

/// @brief The dc motor itself.
struct ObserverSave {
    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> speed;
    std::vector<double> temperature;

    std::size_t decimation;

    ObserverSave()
        : time(),
          current(),
          speed(),
          temperature(),
          decimation(1),
          decimation_cnt(0)
    {
        // Nothing to do.
    }

    void operator()(const State &x, const Time &t) noexcept
    {
        if (++decimation_cnt == decimation) {
            time.emplace_back(t);
            current.emplace_back(x[0]);
            speed.emplace_back(x[1]);
            temperature.emplace_back(x[2]);
            decimation_cnt = 0;
        }
    }

private:
    std::size_t decimation_cnt;
};

/// @brief The dc motor itself.
struct ObserverPrint {
    void operator()(const State &x, const Time &t)
    {
        std::cout << std::fixed << std::setprecision(4) << t << " " << x[0] << " " << x[1] << " " << x[2] << "\n";
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

} // namespace dcmotor_l_temp