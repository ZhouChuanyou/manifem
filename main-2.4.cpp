
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

	auto & hiperbola = RR2.implicit ( x*x + y*y - 1. );

	auto & N = Cell::point();  x ==  0.;  y ==  1.;
	auto & W = Cell::point();  x == -1.;  y ==  0.;
	auto & S = Cell::point();  x ==  0.;  y == -1.;
	auto & E = Cell::point();  x ==  1.;  y ==  0.;

	auto & NW = Mesh::segment ( N, W, 10 );
	auto & WS = Mesh::segment ( W, S, 10 );
	auto & SE = Mesh::segment ( S, E, 10 );
	auto & EN = Mesh::segment ( E, N, 10 );
	auto & circle = Mesh::join ( NW, WS, SE, EN );

	circle.draw_ps ("circle.eps");
  circle.export_msh ("circle.msh");
	
	 cout << "reached end" << endl;
}
