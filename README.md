# PME-DG-Limiters
To compare several limiters 
by using Runge-Kutta Discontinous Galerkin Finite Element Method  
to solve the following Porous Medium Equations
with periodic boundary condition and different initial conditions.
> u_t = (u^m)_xx - c * u^p

## Parameters
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

The parameters are stored in an array named PARA :

```Matlab
PARA = [m, c, p;
		R_left, R_right, dT;
		J, basis_order, 0;
		type_problem, type_limiter, CFL;];
```

## Scripts
### main_test.m
1. Objective:  Test the CFL condition.
2. Output:     CFL_table, which stores the minimum CFL number.
3. Problem:    PME.
4. Initial:    Barenblatt Solution.
5. Solution:   Barenblatt Solution.

### main_ex.m
1. Objective:  Computes the error table.
2. Output:     err_table.
3. Problem:    PME.
4. Initial:    Barenblatt Solution.
5. Solution:   Barenblatt Solution.

### main_app1.m
1. Objective:  Stores the moments of the numerical solution.
2. Output:     Screenshots of the solution.
3. Problem:    PME.
4. Initial:    Two box solutions.
5. Solution:   Unknown.

### main_app2.m
1. Objective:  Stores the moments of the numerical solution.
2. Output:     Screenshots of the solution.
3. Problem:    PME.
4. Initial:    cos(x).
5. Solution:   Unknown.

### main_app3.m
1. Objective:  Stores the moments of the numerical solution.
2. Output:     Screenshots of the solution.
3. Problem:    PME with absorption.
4. Initial:    |Sin(x)| with platform.
5. Solution:   Unknown.

## Functions
### PME.m
This function is the solver of PME and it returns the error.
When type_problem = 0 or 1, it returns the error,
otherwise it returns -1.

There are three global variables which store the information of
limiters.

1. count_mean  where the mean of u is < 0.
2. count_osc   where the oscillation exist.
3. count_pos   where the negative value exist.
### limiter_zq.m
> Numerical Simulation for Porous Medium Equation
> by Local Discontinuous Galerkin Finite Element Method
>
> Auther: Qiang Zhang and Zi-Long Wu