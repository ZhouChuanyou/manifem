
#include "Mesh.h"

using namespace std;

int main () {

   auto environment = Manifold::euclid (2);  // Manifold & environment
   auto & xy = environment.coordinate_system ("Lagrange degree one");
   // FunctionOnMesh::Function & xy
   auto & x = xy[0], & y = xy[1];  // FuncionOnMesh::Function & x, & y
   auto & A = Cell::point ("A");  x ==  0.;  y == 0.;   // Cell & A
   auto & B = Cell::point ("B");  x ==  1.;  y == 0.;
   auto & C = Cell::point ("C");  x ==  0.;  y == 1.;
   auto & M = Cell::point ("M");  x ==  0.5;  y == 0.5;
   auto & AB = Mesh::segment ( A, B, 8 );   // Mesh & AB
   auto & BM = Mesh::segment ( B, M, 5 );
   auto & MC = Mesh::segment ( M, C, 3 );
	 auto & BC = Mesh::join ( BM, MC );
	 auto & CA = Mesh::segment ( C, A, 8 );
   auto & ABC = Mesh::triangle ( AB, BC, CA );   // Mesh & ABC

	 ABC.draw_ps ("tri.eps");
	
	 cout << "reached end" << endl;
}
