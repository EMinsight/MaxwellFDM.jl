# MaxwellFDFD

[![CI](https://github.com/wsshin/MaxwellFDFD.jl/workflows/CI/badge.svg)](https://github.com/wsshin/MaxwellFDFD.jl/actions)
[![Codecov](http://codecov.io/github/wsshin/MaxwellFDFD.jl/coverage.svg?branch=main)](http://codecov.io/github/wsshin/MaxwellFDFD.jl?branch=main)

**MaxwellFDFD** is a high-performance finite-difference frequency-domain (FDFD) solver package
written in Julia.  This package aims to combine the functionalities of
[MaxwellFDFD](https://github.com/wsshin/maxwellfdfd) (written in MATLAB) and
[FD3D](https://github.com/wsshin/fd3d) (written in C using the PETSc library).

Previously, in order to solve large-scale frequency-domain Maxwell's equations on a
finite-difference grid, the users had to go through a tedious procedure as follows:

- Create an "input file" describing the problem using MaxwellFDFD in MATLAB.
- Run FD3D on the input file to solve the problem by iterative methods.
- Load the solution file in MaxwellFDFD for analysis.

Because this procedure involved communication between MaxwellFDFD and FD3D via input and
output files, it was difficult to write scripts that use the solution of the current solve
to create a next problem, which is the capability needed in, e.g., the inverse design
procedure.  Also, because users often do not have a MATLAB license on the computation
cluster, the above procedure typically involves uploading input/downloading output files
between a local computer with a MATLAB license and computation cluster, thereby making the
situation even worse.


MaxwellFDFD aims to avoid these complications by implementing an FDFD solver in Julia.  Julia
can easily replace any programs written in MATLAB, so it is straightforward to implement the
capabilities of MaxwellFDFD in MaxwellFDFD.  Furthermore, using the [Julia wrapper of the
PETSc library](https://github.com/JuliaParallel/PETSc.jl), we can implement the capabilities
of FD3D, i.e., iterative solution algorithms for distributed—memory computation clusters,
directly within MaxwellFDFD.  This means that you don't need to deal with multiple programs
to launch a next solve based on the current solution.

In addition, MaxwellFDFD aims to implement capabilities lacking in MaxwellFDFD and FD3D, such
as

- subpixel smoothing of material parameters
- symmetry boundary on the positive end of the computation domain
- more general TF/SF method

MaxwellFDFD is still a project under development.

<!---Mention Ian Williamson's repository.--->
