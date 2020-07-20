# step-amr: Parallel Mesh Refinement with p4est

This example is a new example showing how to construct and refine meshes in
parallel. The parallel::distributed::Triangulation internally uses p4est for
the parallel partitioning. The tutorial is a simplification of [deal.II
step-40](https://www.dealii.org/current/doxygen/deal.II/step_40.html).

## What it does

- Run in parallel using MPI
- Represent mesh in parallel using p4est
- Interpolate a function into a Finite Element space
- Refine the mesh based on a simple error estimator
- Output the solution in parallel

## Exercises

1. Open up and visualize the solution using the "warp by scalar" filter in
   ParaView. Why are there small gaps where adaptive refinement is being done?
   Hint: switch to linear finite elements (degree=1). Do they disappear?

2. What needs to be changed to run the program in 3d? Visualize the solution
   in ParaView!

3. Change the geometry: maybe a unit circle or a torus? see the [GridGenerator
   namespace](https://www.dealii.org/9.2.0/doxygen/deal.II/namespaceGridGenerator.html)
   for inspiration. Further,
   [step-49](https://www.dealii.org/current/doxygen/deal.II/step_49.html)
   documents advanced way to create and manipulate meshes.

