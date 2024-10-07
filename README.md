# FEM_triangles_order1

Educational python code that have to be filled. It uses the finit element method to solve the heat equation or the helmholtz one.

# Initial commit

```
git clone [URL]
```
If you want to create meshes, you can use the file [Meshes/msh_creation.py](Meshes/msh_creation.py). You will need the gmsh library, with the python binding.

# Files to change

- [solve_helmholtz.py](solve_helmholtz.py) : define the manufactured solution g(x,y) and the corresponding source term f(x,y)
- [FEMlib/errors.py](FEMlib/errors.py) : add the code corresponding to the computation of the error in L2 norm

