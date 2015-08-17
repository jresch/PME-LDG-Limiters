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
Initial:    Barenblatt Solution.
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
This function is the solver of PME and it returns the error.
When type_problem = 0 or 1, it returns the error,
otherwise it returns -1.

There are three global variables which tracks the limiters.

- track_mean  where the mean of u is < 0.
- track_osc   where the oscillation exist.
- track_pos   where the negative value exist.

### limiter_zq.m
> Numerical Simulation for Porous Medium Equation
> by Local Discontinuous Galerkin Finite Element Method
>
> Auther: Qiang Zhang and Zi-Long Wu