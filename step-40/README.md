# step-40: Solve Laplace in parallel

This example is based on [deal.II step-40](https://www.dealii.org/current/doxygen/deal.II/step_40.html). Please go there for additional information.

## What it does

- Run in parallel using MPI
- Represent mesh in parallel using p4est
- Solve a Laplace problem using algebraic multigrid
- Using PETSc or Trilinos for linear algebra
- Refine the mesh adaptively based on a simple error estimator

## Exercises

1. Add text output that shows the number of locally owned cells and locally owned DoFs on each processor.

