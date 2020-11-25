# How ricpad sets up precision

There are three important variables that determine how accurate a ricpad calculation is:
* digits: sets the amount of digits for the numerical representations.
* nr_tolerance: sets the difference in the eigenvalue between two Newton-Raphson steps.
* nr_diff_epsilon: sets the step size for the Newton-Raphson method

First, it checks if target_digits has been set. If so, nr_tolerance is set to 10^(-target_digits - 20), nr_step_size is set to 10^(-(nr_step_size + 100)), and digits is set to nr_step_size + 200.

If target_digits hasn't been defined, ricpad checks if digits has been defined; if so, nr_step_size is set to 10^(-digits+100), and nr_tolerance is set to 10^(-digits - 120).
