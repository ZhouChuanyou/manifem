
#include "Mesh.h"
#include "math.h"

using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	auto RR2 = Manifold::euclid (2);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'R3' :
	auto & xy = RR2.coordinate_system ("Lagrange degree one");

	// we can extract components of xy using the [] operator :
	FunctionOnMesh::Function & x = xy[0], & y = xy[1];

	auto & hiperbola = RR2.implicit ( x*y - 1. );

	auto & A = Cell::point();  x ==  0.5;   y ==  2.;
	auto & B = Cell::point();  x ==  3;     y ==  0.333333333333;

	auto & arc_of_hiperbola = Mesh::segment ( A, B, 15 );

	// arc_of_hiperbola.draw_ps ("hiperbola.eps");
  arc_of_hiperbola.export_msh ("hiperbola.msh");
	
	 cout << "reached end" << endl;
}
