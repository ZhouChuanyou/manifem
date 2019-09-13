
#include "Mesh.h"
#include "math.h"

using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	auto RR3 = Manifold::euclid (3);
	
	// xyz is a map defined on our future mesh with values in 'environment' :
	auto & xyz = RR3.coordinate_system ("Lagrange degree one");
	// we can extract components of xyz using the [] operator :
	FunctionOnMesh::Function & x = xyz[0], & y = xyz[1], & z = xyz[2];

	auto & sphere = RR3.implicit ( x*x + y*y + z*z - 1. );

	// let's mesh half of a sphere
	auto & E  = Cell::point();  x ==  1.;   y ==  0.;   z == 0.;
	auto & N  = Cell::point();  x ==  0.;   y ==  1.;   z == 0.;
	auto & W  = Cell::point();  x == -1.;   y ==  0.;   z == 0.;
	auto & S  = Cell::point();  x ==  0.;   y == -1.;   z == 0.;
	auto & up = Cell::point();  x ==  0.;   y ==  0.;   z == 1.;
	int n = 15;
	auto & EN = Mesh::segment ( E, N, n );
	auto & NW = Mesh::segment ( N, W, n );
	auto & WS = Mesh::segment ( W, S, n );
	auto & SE = Mesh::segment ( S, E, n );
	auto & upE = Mesh::segment ( up, E, n );
	auto & upN = Mesh::segment ( up, N, n );
	auto & upW = Mesh::segment ( up, W, n );
	auto & upS = Mesh::segment ( up, S, n );
	
	// now four triangles
	auto & ENup = Mesh::triangle ( EN, upN.reverse(), upE );
	auto & NWup = Mesh::triangle ( NW, upW.reverse(), upN );
	auto & WSup = Mesh::triangle ( WS, upS.reverse(), upW );
	auto & SEup = Mesh::triangle ( SE, upE.reverse(), upS );

	// and finally join the triangles :
	auto & hemisphere = Mesh::join ( ENup, NWup, WSup, SEup );
	
  hemisphere.export_msh ("hemisphere.msh");
	
	cout << "reached end" << endl;
}
