# deal.II - p4est Tutorial

This is the example repository for the 2020 p4est&amp;deal.II tutorial. The
goal is to teach basics about the [deal.II library](https://dealii.org) and
how it interfaces with [p4est](http://p4est.org). For more information, please
see the tutorial website at http://www.p4est.org/school.html.


## Prerequisites

To run the examples you will need:
- deal.II version 9.2.0 or newer
- installed with dependencies p4est, Trilinos or PETSc
- ParaView to visualize solutions

Ways to run the tutorials:
1. VirtualBox image at https://www.math.clemson.edu/~heister/dealvm/
2. Docker images at https://github.com/dealii/docker-files (if you are experienced with using Docker)
3. A manual installation, see [deal.II download page](https://www.dealii.org/download.html)

## Additional Resources

Some useful resources to study:
- https://www.dealii.org/ - The main website of deal.II
- [deal.II tutorials](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html) - Tutorials and Code documentation
- The original [deal.II and p4est paper](https://dl.acm.org/doi/10.1145/2049673.2049678) ([preprint](http://www.math.clemson.edu/~heister/preprints/BangerthBursteddeHeisterKronbichler_distributed.pdf))
- [Wolfgang's video lectures](https://www.math.colostate.edu/~bangerth/videos.html)


## The tutorials

1. [``step-2/`` Meshes and Finite Element spaces](./step-2/) - slightly modified [deal.II step-2](https://www.dealii.org/current/doxygen/deal.II/step_2.html)
2. [``step-amr/`` Parallel Mesh Refinement with p4est](./step-amr/) - a new program
3. [``step-40/`` Solve Laplace in parallel](./step-40/) - slightly modified [deal.II step-40](https://www.dealii.org/current/doxygen/deal.II/step_40.html)
4. Take a look at [deal.II step-50](https://www.dealii.org/current/doxygen/deal.II/step_50.html): large-scale matrix-free geometric multigrid.
