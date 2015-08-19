# PME-DG-Limiters
To compare several limiters
by using **Runge-Kutta Discontinous Galerkin Finite Element Method**
to solve the following **Porous Medium Equation**
with **periodic boundary condition** and different initial conditions.
> u_t = (u^m)_xx - c * u^p

## Flow Chart

```
main_xx.m
  - Define the paths of data and report.
  - Define the parameters of problem (equation, boundary, initial).
  - Define the parameters of method (cells, basis, limiter).
  - Call PME.m to solve the problem.
      - Extract parameters.
      - Call Quadrature_Set.m to set up quadrature.
      - Call Mesh_Set.m to generate mesh.
      - Define space and time step
        according to the pre-stored CFL condition number.
      - Call Basis_Set.m to set up basis.
      - Initialize the variables to track limiter.
      - Set up the initial condition.
      - Execute the for loop to solve the problem
        > Using 3rd order Runge-Kutta time marching.
        > For applications, the moments of solutions will be saved.
        - Call L_pme.m to march.
            - Call H_pme.m to calculate the flux.
        - Call Limiter.m to adjust the solution.
          - Call the corresponding limiter.
      - If using limiter, save the trackers and other variables.
```

## Data Structure
### Parameters
Parameters for solving the problem.

```
Parameters				| Description
----------------------- | -----------
m, c, p					| Parameters of equation.
R_left, R_right         | Region is [R_left, R_right].
dT                      | Time duration.
uI                      | Initial value of solution.
J                       | Number of cells.
basis_order             | Basis order.
type_problem            | 0: CFL test (call main_test).
                        | 1: Experiment (call main_ex).
                        | 2: Application (call main_app).
type_limiter            | Different limiters.
CFL                     | In CFL test, CFL is given outside PME.m.
mu = 1                  | Parameter of the threshold in the limiter
                        | who attenuates the oscillation.
num_fig = 200           | Number of figures of the solution.
```

### Variables
Variables during solving the PME problem.

```
Variables				| Description
----------------------- | -----------
PARA                    | PARA = [m, c, p;
                        |         R_left, R_right, dT;
                        |         J, basis_order, 0;
                        |         type_problem, type_limiter, CFL;];
                        | Scripts use PARA to transmit parameters.
points, weights         | size(points) = (7, 1)
                        | size(weights) = (1, 7)
                        | The gauss points in [-1, 1] and the corresponding weights.
                        | (Using five points Gauss–Legendre quadrature.)
grid                    | size(grid) = (2, J)
                        | grid(1, j) stores the center of the j-th cell.
                        | grid(2, j) stores half of the length of the j-th cell.
X                       | size(X) = (7, J)
                        | X(:, j) stores the guass points in the j-th cell.
dx                      | Space step.
dt                      | Time step.
loops                   | Number of loops.
psi                     | size(psi) = (7, basis_order + 1)
                        | psi(:, i + 1) stores the values of the i-th degree Legendre polynomial
                        | on the gauss points.
psi_z                   | size(psi_z) = (7, basis_order + 1)
                        | psi_z(:, i + 1) stores the values of the drivative
                        | of the i-th degree Legendre polynomial on the gauss points.
u                       | size(u) = (7, J)
                        | u = u(X), X(:, j) stores the values of solution on the guass points in the j-th cell.
u_coord                 | size(u_coord) = (basis_order + 1, J)
                        | u_coord(i + 1, j) stores the weights of the i-th degree Legendre polynomial.
track_mean              | size(track_mean) = (J, loops), it's a sparse matrix.
                        | track_mean(j, loop) = cell average,
                        | if the cell average of the j-th cell is negative, otherwise it equals 0.
track_osc               | size(track_osc) = (J, loops), it's a sparse matrix.
                        | track_osc(j, loop) ~= 0, if there's oscillation in the j-th cell.
track_pos               | size(track_pos) = (J, loops), it's a sparse matrix.
                        | track_pos(j, loop) ~= 0, if there's negative value in the j-th cell.
```

## Scripts
### main_test.m

```
Objective:  Test the CFL condition.
Output:     CFL_table, which stores the minimum CFL number.
Problem:    PME.
Initial:    Barenblatt Solutison.
Solution:   Barenblatt Solution.
```

### main_ex.m

```
Objective:  Computes the error table.
Output:     err_table.
Problem:    PME.
Initial:    Barenblatt Solution.
Solution:   Barenblatt Solution.
```

When writing new limiters, we can use this script to test
if the limiter works and how the limiter works.

### main_app1.m

```
Objective:  Stores the moments of the numerical solution.
Output:     Screenshots of the solution.
Problem:    PME.
Initial:    Two box solutions.
Solution:   Unknown.
```

### main_app2.m

```
Objective:  Stores the moments of the numerical solution.
Output:     Screenshots of the solution.
Problem:    PME.
Initial:    cos(x).
Solution:   Unknown.
```

### main_app3.m

```
Objective:  Stores the moments of the numerical solution.
Output:     Screenshots of the solution.
Problem:    PME with absorption.
Initial:    |Sin(x)| with platform.
Solution:   Unknown.
```

### Data_Process.m

```
Objective:  Plots the trackers to learn where and how the limiter works.
Output:     The locations of cells where the limiter is used.
            The min and mean of the negative cell averages.
```

## Functions
### PME.m

```Matlab
function [err] = PME(PARA, uI, path_data, path_report)
```
This function is the solver of PME and it returns the error.
When `type_problem` = 0 or 1, it returns the error,
otherwise it returns -1.

This function also saves the moments of the solution to
observe the movement of the solution
and it saves three global variables which tracks the limiters.

- `track_mean`  records the cells where the cell average is negative.
- `track_osc`   records the cells where the oscillation exist.
- `track_pos`   records the cells where the negative value exist.

### L_pme.m

```Matlab
function [ut_coord] = L_pme(u_coord, loop, type_limiter)
```

This function is to solve the L(u)
where u_t = L(u) is the first order ODE system.

If calling limiter_yy.m, this function will call twice H_pme.m
and use the combination of the flux to adjust the cell average
and make sure the cell average is non-negative.
(Sometimes there are some extreme small negative values due to the machine error.
Those values will be replaced with zero and it won't affect the accuracy.)

### H_pme.m

```Matlab
function [flux_ur, q] = H_pme(u, loop, psi, psi_z)
```

This function is to calculate the flux according to the given basis.

### Limiter.m

```Matlab
function [u_coord] = Limiter(u_coord, loop, type_limiter)
```

This function calls the corresponding limiter to adjust the solution.

### limiter_zq.m

```Matlab
function [u_coord] = limiter_zq(u_coord, loop)
```

This limiter is designed according to the following article.
> Numerical Simulation for Porous Medium Equation
> by Local Discontinuous Galerkin Finite Element Method
>
> Auther: Qiang Zhang and Zi-Long Wu

### limiter_yy.m

```Matlab
function [u_coord] = limiter_yy(u_coord, loop)
```

This limiter is designed according to the following article.
> Positivity-preserving high-order local discontinuous Galerkin method
> for parabolic equations with blow-up solutions
>
> Auther: LiGuo and YangYang

This function uses zero to replace the extreme small negative cell averages
and adjust those points with negative values.
The codes of adjusting negative cell average is in L_pme.m.

### limiter_pos.m

```Matlab
function [u_coord] = limiter_pos(u_coord, loop)
```

This limiter uses 2 degree polynomial to attenuate the oscillation
and uses the same process as limiter_zq to keep solution positive.

### Quadrature_Set.m

```Matlab
function [points, weights] = Quadrature_Set()
```

This function generates the gauss points and the corresponding weights
for the quadrature.

`points` is column vector and `weights` is row vector.

```
\int\limits_{-1}^{1}f(x)dx = weights * f(points)
```

### Mesh_Set.m

```Matlab
function [grid,X] = Mesh_Set(points,R_left,R_right,J)
```

This function generates the mesh. (mesh is a reserved word in Matlab.)

`grid` stores the center and half length of each cell.
`X` is the mesh.

### Basis_Set.m

```Matlab
function [psi,psi_z] = Basis_Set(points,basis_order)
```

This function returns the values of basis
on the gauss points according to the order of basis.

### BarenblattSolution.m

```Matlab
function [y] = BarenblattSolution(x, t, m)
```

This function is the famous Barenblatt Solution
which is the exact solution of Porous Medium Equation.

### minmod.m

```Matlab
function [a, flag] = minmod(a1, a2, a3, threshold)
```

This minmod function is used in limiters. The principle is:

```
minmod(a1, a2, a3, threshold)
    = a1                            if abs(a1) <= threshold
      s * min(abs([a1, a2, a3]))    if s = sign(ai), i = 1, 2, 3
      0                             otherwise
```