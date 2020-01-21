# manifem
ManiFEM is a software for solving partial differential equations through the finite element method.
The name comes from "finite elements on manifolds". 
ManiFEM has been designed with the goal of coping with very general meshes,
in particular meshes on Riemannian manifolds, even manifolds which cannot be embedded in R^3, like the torus R^2/Z^2.
Also, ManiFEM has been written with the goal of being conceptually clear and easy to read.
We hope it will prove particularly useful for people who want a fine control over the mesh, 
e.g. for implementing their own meshing or remeshing algorithms.

ManiFEM is written in C++ and uses the Standard Template Library (STL). 
It uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for storing matrices and solving systems of linear equations.

ManiFEM is just a collection of C++ classes. It has no user-friendly interface nor graphic capabilities. The user shold have some understanding of programming and of C++. However, maniFEM can be used at a basic level by people with no deep knowledge of C++. 

ManiFEM is being developped by Cristian Barbarosie, Sérgio Lopes and Anca-Maria Toader;
see its [homepage](https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/).

To learn maniFEM, you should begin by reading the [manual](https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/manual-manifem.pdf).

To use maniFEM, you need a recent version of the `g++` compiler and the `make` utility. Under linux it should be easy to install them. It is not that easy to install and use them under windows, but it is certainly possible, for instance by using [cygwin](https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/cygwin.org). You may also want to use [gmsh](http://gmsh.info/) for visualization purposes. You should copy all `.cpp` and `.h` files to some directory in your computer, together with the `Makefile`.
Then try `make run-1.1` for the example in paragraph 1.1 of the manual, `make run-2.2` for the example in paragraph 2.2, and so on.

