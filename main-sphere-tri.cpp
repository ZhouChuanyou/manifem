
#include "Mesh.h"
#include "math.h"

using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	auto RR3 = Manifold::euclid (3);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'environment' :
	auto & xyz = RR3.coordinate_system ("Lagrange degree one");
	xyz.components[0]->name = "x";
	xyz.components[1]->name = "y";
	xyz.components[2]->name = "z";

	// we can extract components of xyz using the [] operator :
	FunctionOnMesh::Function & x = xyz[0], & y = xyz[1], & z = xyz[2];

	auto & sphere = R3.implicit ( x*x + y*y + z*z + 1.5*x*y*z - 1. );

	// let's mesh half of a sphere
	auto & A = Cell::point();  x ==  1.;   y ==  0.;   z == 0.;
	auto & B = Cell::point();  x ==  0.;   y ==  1.;   z == 0.;
	auto & C = Cell::point();  x == -1.;   y ==  0.;   z == 0.;
	auto & D = Cell::point();  x ==  0.;   y == -1.;   z == 0.;
	auto & N = Cell::point();  x ==  0.;   y ==  0.;   z == 1.;
	int n = 25;
	auto & AB = Mesh::segment ( A, B, n );
	auto & BC = Mesh::segment ( B, C, n );
	auto & CD = Mesh::segment ( C, D, n );
	auto & DA = Mesh::segment ( D, A, n );
	auto & NA = Mesh::segment ( N, A, n );
	auto & NB = Mesh::segment ( N, B, n );
	auto & NC = Mesh::segment ( N, C, n );
	auto & ND = Mesh::segment ( N, D, n );
	
	// now four triangles
	auto & ABN = Mesh::triangle ( AB, NB.reverse(), NA );
	auto & BCN = Mesh::triangle ( BC, NC.reverse(), NB );
	auto & CDN = Mesh::triangle ( CD, ND.reverse(), NC );
	auto & DAN = Mesh::triangle ( DA, NA.reverse(), ND );

	// and finally join the triangles :
	auto & hemisphere = Mesh::join ( ABN, BCN, CDN, DAN );
	
  hemisphere.export_msh ("hemisphere.msh");
	
	 cout << "reached end" << endl;
}
