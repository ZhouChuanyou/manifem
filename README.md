# manifem
ManiFEM is a software for solving partial differential equations through the finite element method.
The name comes from "finite elements on manifolds". 
ManiFEM has been designed with the goal of coping with very general meshes,
in particular meshes on Riemannian manifolds, even manifolds which cannot be embedded in R^3, like the torus R^2/Z^2.
Also, ManiFEM has been written with the goal of being conceptually clear and easy to read.

ManiFEM is written in C++ and uses the Standard Template Library (STL). 
It uses Eigen for storing matrices and solving systems of linear equations, eigen.tuxfamily.org/index.php?title=Main_Page.

To learn ManiFEM, you should probably begin by reading the document manual.pdf.

To use ManiFEM, copy all .cpp and .h files to some directory in your computer, together with Makefile.
Then try 'make run-1.1' for the example in paragraph 1.1 of the manual, 'make run-2.2' for the example in paragraph 2.2, and so on.

