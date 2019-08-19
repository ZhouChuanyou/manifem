
#include "Mesh.h"

using namespace std;


int main () {

	// we choose our (geometric) space dimension :
	auto & environment = Manifold::euclid (3);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'environment' :
	auto & xyz = environment.coordinate_system ("Lagrange degree one");

	// we can extract components of xyz using the [] operator :
	auto & x = xyz[0], & y = xyz[1], & z = xyz[2];

	// Let's build a rectangular mesh. First, the four corners :
	auto & SW = Cell::point ("SW");  x ==  0.;  y ==  3.;
	auto & SE = Cell::point ("SE");  x ==  2.;  y ==  3.;
	auto & NE = Cell::point ("NE");  x ==  2.;  y ==  5.;
	auto & NW = Cell::point ("NW");  x ==  0.;  y ==  5.;
	// Now build the four sides of the rectangle :
	auto & south = Mesh::segment ( SW, SE, 2 );
	auto & east  = Mesh::segment ( SE, NE, 2 );
	auto & north = Mesh::segment ( NE, NW, 2 );
	auto & west  = Mesh::segment ( NW, SW, 2 );
	// and now the rectangle :
	// auto & rect_mesh = Mesh::rectangle ( south, east, north, west );
	auto & rect_mesh = Mesh::rectangle ( south, east, north, west, tag::with_triangles );
	
	Cell & P = Cell::point ("P"); z == 6.;
	Cell & R = Cell::point ("R"); z == 8.;
	
	Mesh & PR = Mesh::segment ( P, R, 2 );
	
	Mesh & product_mesh = Mesh::cartesian_product ( rect_mesh, PR );

	product_mesh.export_msh ("cartesian-tri.msh");

} // end of main
