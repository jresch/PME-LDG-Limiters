# PME-DG-Limiters
To compare several limiters
by using Runge-Kutta Discontinous Galerkin Finite Element Method
to solve the following Porous Medium Equations
with periodic boundary condition and different initial conditions.
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
        - Call Limiter.m to adjust the solution.
          - Call the corresponding limiter.
      - If using limiter, save the trackers and other variables.
```

## Parameters

```
Variables				| Description
----------------------- | -----------
m, c, p					| Equation
R_left, R_right, dT     | Boundary
uI                      | The initial value
J                       | Number of cells
basis_order             | Basis
type_problem            | 0: CFL test (call main_test)
                        | 1: experiment (call main_ex)
                        | 2: application (call main_app)
type_limiter            | Different limiters
CFL                     | In CFL test, CFL is given outside PME.m.
```
The parameters are stored in an array named PARA :

```Matlab
PARA = [m, c, p;
		R_left, R_right, dT;
		J, basis_order, 0;
		type_problem, type_limiter, CFL;];
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

## Functions
### PME.m

```Matlab
function [err] = PME(PARA, uI, path_data, path_report)
```
This function is the solver of PME and it returns the error.
When type_problem = 0 or 1, it returns the error,
otherwise it returns -1.

This function also saves the moments of the solution to
observe the movement of the solution
and it saves three global variables which tracks the limiters.

- track_mean  records the cells where the cell average is negative.
- track_osc   records the cells where the oscillation exist.
- track_pos   records the cells where the negative value exist.

### L_pme.m

```Matlab
function [u_coord] = L_pme(u_coord, loop)
```

This function is to solve the L(u)
where u_t = L(u) is the first order ODE system.

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

### Quadrature_Set.m

```Matlab
function [points, weights] = Quadrature_Set()
```

This function generates the gauss points and the corresponding weights
for the quadrature.

*points* is column vector and *weights* is row vector.

```
\int\limits_{-1}^{1}f(x)dx = weights * f(points)
```

### Mesh_Set.m

```Matlab
function [grid,X] = Mesh_Set(points,R_left,R_right,J)
```

This function generates the mesh. (mesh is a reserved word in Matlab.)

*grid* stores the center and half length of each cell.
*X* is the mesh.

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