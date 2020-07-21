# step-2: Create and adapt a mesh and Finite Element spaces

This example is based on [deal.II
step-2](https://www.dealii.org/current/doxygen/deal.II/step_2.html). Please
see step-1 and step-2 there for more information.

## What it does

- This is a sequential code (no MPI)
- We adaptively refine a mesh by iterating over cells
- Distribute degrees of freedom (DoFs) of a Finite Element space
- Generate and print sparsity patterns for a sparse matrix
- Renumber degrees of freedom

## Exercises

1. Get the source:
```
    git clone https://github.com/tjhei/dealii-p4est-tutorial
```
2. Run the program:
```
    cmake .
    make
    ./step-2
```
3. Look at the sparsity pattern and see how it changes when you change the
   finite element degree (set in the constructor of FE_Q).

4. Do the visualization of the DoF and support points under [step-2 Possible
   extensions](https://www.dealii.org/9.2.0/doxygen/deal.II/step_2.html#Possibleextensions). Maybe with a simpler mesh...
