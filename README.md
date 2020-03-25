# maniFEM
ManiFEM is a C++ library for solving partial differential equations through the finite element method.
The name comes from "finite elements on manifolds". 
ManiFEM has been designed with the goal of coping with very general meshes,
in particular meshes on Riemannian manifolds, even manifolds which cannot be embedded in R^3, like the torus R^2/Z^2.
Also, maniFEM has been written with the goal of being conceptually clear and easy to read.
We hope it will prove particularly useful for people who want a fine control over the mesh, 
e.g. for implementing their own meshing or remeshing algorithms.

ManiFEM uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for storing matrices and for
solving systems of linear equations.

ManiFEM is just a collection of C++ classes. It has no user-friendly interface nor graphic capabilities. 
The user shold have some understanding of programming and of C++. 
However, maniFEM can be used at a basic level by people with no deep knowledge of C++.

In its current version (20.03), maniFEM works well for mesh generation. 
Quotient manifolds and anisotropic Riemann metrics are not yet implemented. 
Variational formulations and finite elements are not yet implemented.
To check which version of maniFEM is installed in your computer, see at the beginning of the file `maniFEM.h`.

A component of maniFEM, [MetricTree](https://github.com/cristian-barbarosie/MetricTree), can be used independently.

ManiFEM is being developped by [Cristian Barbarosie](cristian.barbarosie@gmail.com), Sérgio Lopes and Anca-Maria Toader;
see its [homepage](https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/).

To learn maniFEM, you should read the 
[manual](https://webpages.ciencias.ulisboa.pt/~cabarbarosie/manifem/manual-manifem.pdf).

To use maniFEM, you should copy all files under [src/](https://github.com/cristian-barbarosie/manifem/blob/master/src/) to some directory in your computer.
Then try `make run-1.1` for the example in paragraph 1.1 of the manual, `make run-2.2` for the example in paragraph 2.2, 
and so on.
You will need a recent C++ compiler (we use `g++`) and the `make` utility. 
Under linux it should be easy to install them. 
It is not that easy to install and use them under Windows, but it is certainly possible, for instance by using 
[cygwin](https://cygwin.org). 
You may also want to use [gmsh](http://gmsh.info/) for visualization purposes. 

This work is supported by National Funding from FCT - Fundação para a Ciência e a Tecnologia (Portugal), 
through Faculdade de Ciências da Universidade de Lisboa and 
Centro de Matemática, Aplicações Fundamentais e Investigação Operacional, project UID/MAT/04561/2020.

Copyright 2019, 2020 Cristian Barbarosie cristian.barbarosie@gmail.com

ManiFEM is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ManiFEM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

Full text of the GNU Lesser General Public License is available in files [COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER).
It can also be found at <https://www.gnu.org/licenses/>.
