# How to use ricpad

The program ricpad takes an input file that needs to have the following entries:

```
pot <potential> # Expression for the potential. For example, x^4
var <variable>  # Name of the variable used to define the potential.
D <Dini> <Dmax> # Value of the initial and final sizes of the Hankel matrices
d <d>           # Value of the RPM displacement parameter, d.
s <s>           # Value of the parity parameter for even problems, or for the angular momentum for radial problems
E0 <E0>         # Initial value of the energy for the NR method
```

Additionally, it can define the following options:

```
use_complex <0 or 1>          # Use mpc_complex instead of mpfr_float. 1 (True) by default.
problem_type <even or radial> # Use the RPM to solve even or radial problems. even by default
target_digits <target_digits> # Amount of digits we want the result to be accurate. 40 by default.
nr_step_size <nr_step_size>   # Size of the h value in the NR method to compute numerical differentiation. See following section for description.
nr_tolerance <nr_tolerance>   # Difference between steps of the NR method to consider convergence reached.
E0I <E0I>                     # Value of the imaginary part of the initial energy. 0 by default.
```
To call `ricpad`, run:
```
ricpad [input_file]
```

See the examples folder for some examples.

# How ricpad sets up precision

There are three important variables that determine how accurate a ricpad calculation is:
* `digits`: sets the amount of digits for the numerical representations.
* `nr_tolerance`: sets the difference in the eigenvalue between two Newton-Raphson steps.
* `nr_diff_epsilon`: sets the step size for the Newton-Raphson method

First, it checks if target_digits has been set. If so, `nr_tolerance` is set to `10^(-target_digits - 20)`, `nr_step_size` is set to `10^(-(nr_step_size + 100))`, and digits is set to `nr_step_size + 200`.

If target_digits hasn't been defined, ricpad checks if digits has been defined; if so, `nr_step_size` is set to `10^(-digits+100)`, and `nr_tolerance` is set to `10^(-digits - 120)`.
