---
layout: default
---

# NumInt

[![Ubuntu](https://github.com/Galfurian/numint/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/Galfurian/numint/actions/workflows/ubuntu.yml)
[![Windows](https://github.com/Galfurian/numint/actions/workflows/windows.yml/badge.svg)](https://github.com/Galfurian/numint/actions/workflows/windows.yml)
[![MacOS](https://github.com/Galfurian/numint/actions/workflows/macos.yml/badge.svg)](https://github.com/Galfurian/numint/actions/workflows/macos.yml)
[![Documentation](https://github.com/Galfurian/numint/actions/workflows/documentation.yml/badge.svg)](https://github.com/Galfurian/numint/actions/workflows/documentation.yml)

A high-performance numerical integration library inspired by
[odeint-v2](https://github.com/headmyshoulder/odeint-v2). This library provides
tools for solving ordinary differential equations (ODEs) using various
integration methods, including adaptive and fixed-step approaches.

## Features

- **Adaptive and Fixed-Step Integration**:
  - Adaptive step-size control for efficient and accurate simulations.
  - Fixed-step integration for simple and predictable simulations.
- **Integration Methods**:
  - Euler Method
  - Improved Euler Method (Heun's Method)
  - Runge-Kutta 4th Order (RK4)
- **Customizability**:
  - Support for user-defined termination conditions.
  - Decimation for efficient observation.
- **Error Control**:
  - Absolute, relative, and mixed truncation error handling.

## Getting Started

### Installation

Clone the repository and include the `numint` folder in your project's include path.

```bash
git clone https://github.com/yourusername/numerical-integration-library.git
```

### Prerequisites

- C++17 or later.
- (Automatically fetched) [gpcpp](https://github.com/yourusername/gpcpp) for visualization.

## Usage

### Example: Fixed-Step Integration

```cpp
#include <numint/solver.hpp>
#include <numint/stepper/stepper_euler.hpp>
#include <iostream>
#include <array>

/// Define the state.
using State = std::array<double, 2>;

/// Define the system model.
class Model {
public:
    void operator()(const State &state, State &dxdt, double) const
    {
        dxdt[0] = state[1];
        dxdt[1] = -state[0];
    }
};

int main()
{
    // Initial state: [position, velocity].
    State state = { 1.0, 0.0 };
    // Instantiate the model.
    Model model;
    // Define and instantiate the stepper.
    numint::stepper_euler<State, double> solver;
    // Simulation parameters.
    const double start_time = 0.0;
    const double end_time   = 10.0;
    const double time_step  = 0.01;
    // Integrate using fixed stepper.
    numint::integrate_fixed(
        solver,
        [](const State &x, double t) {
            std::cout << "Time: " << t << ", State: [" << x[0] << ", " << x[1] << "]\n";
        },
        model,
        state,
        start_time,
        end_time,
        time_step);
    return 0;
}
```

if you run it, it will produce:

```bash
-> % ./numint_dcmotor
Time: 0, State: [1, 0]
Time: 0, State: [1, -0.01]
Time: 0.01, State: [0.9999, -0.02]
...
Time: 9.99, State: [-0.88228, 0.571618]
Time: 10, State: [-0.876564, 0.580441]
Number of integration points :1001
```

### Example: Adaptive-Step Integration

```cpp
#include <numint/solver.hpp>
#include <numint/stepper/stepper_adaptive.hpp>
#include <numint/stepper/stepper_rk4.hpp>
#include <iostream>
#include <array>

/// Define the state.
using State = std::array<double, 2>;

/// Define the system model.
class Model {
public:
    void operator()(const State &state, State &dxdt, double) const
    {
        dxdt[0] = state[1];
        dxdt[1] = -state[0];
    }
};

int main()
{
    // Initial state: [position, velocity].
    State state = { 1.0, 0.0 };
    // Instantiate the model.
    Model model;
    // Define and instantiate the adaptive stepper.
    numint::stepper_adaptive<numint::stepper_rk4<State, double>> solver;
    // Set adaptive stepper parameters.
    solver.set_tollerance(1e-5);
    solver.set_min_delta(1e-9);
    solver.set_max_delta(0.1);
    // Simulation parameters.
    const double start_time        = 0.0;
    const double end_time          = 10.0;
    const double initial_time_step = 0.01;
    // Integrate using fixed stepper.
    auto points = numint::integrate_adaptive(
        solver,
        [](const State &x, double t) {
            std::cout << "Time: " << t << ", State: [" << x[0] << ", " << x[1] << "]\n";
        },
        model,
        state,
        start_time,
        end_time,
        initial_time_step);
    std::cout << "Number of integration points :" << points << "\n";
    return 0;
}
```

if you run it, it will produce:

```bash
-> % ./numint_dcmotor
Time: 0, State: [0.99995, -0.00999983]
Time: 0.01, State: [0.999608, -0.0279963]
...
Time: 9.91872, State: [-0.839072, 0.544021]
Time: 10, State: [-0.839072, 0.544021]
Number of integration points :103
```

### Visualization

If `ENABLE_PLOT` is defined and [gpcpp](https://github.com/yourusername/gpcpp) is available, results can be visualized using Gnuplot.

## Documentation

### Core Functions

#### `integrate_fixed`

Integrates a system over a fixed time step.

```cpp
int integrate_fixed(Stepper &stepper, Observer &&observer, System &&system, 
                    Stepper::state_type &state, Stepper::time_type start_time,
                    Stepper::time_type end_time, Stepper::time_type time_delta);
```

#### `integrate_adaptive`

Integrates a system using an adaptive stepper.

```cpp
int integrate_adaptive(Stepper &stepper, Observer &&observer, System &&system, 
                       Stepper::state_type &state, Stepper::time_type start_time,
                       Stepper::time_type end_time, Stepper::time_type time_delta);
```

### Available Steppers

The basic steppers:

- `stepper_euler`: Implements Euler's method.
- `stepper_improved_euler`: Implements Heun's method.
- `stepper_midpoint`: Implements the Midpoint method.
- `stepper_rk4`: Implements Runge-Kutta 4th Order.
- `stepper_simpsons`: Implements Simpson's rule for integration.
- `stepper_trapezoidal`: Implements the Trapezoidal rule for integration.

and the adaptive one, which wraps one of the previous steppers:

- `stepper_adaptive`: Dynamically adjusts step size for accuracy and efficiency.

## Contributing

Contributions are welcome! Please submit issues or pull requests to improve the library.

## License

This project is licensed under the MIT License - see the `LICENSE.md` file for details.

---

Enjoy efficient numerical integration with `numint`!