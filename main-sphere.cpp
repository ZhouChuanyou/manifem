
#include "Mesh.h"
#include "math.h"

using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	auto R3 = Manifold::euclid (3);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'environment' :
	auto & xyz = R3.coordinate_system ("Lagrange degree one");
	xyz.components[0]->name = "x";
	xyz.components[1]->name = "y";
	xyz.components[2]->name = "z";

	// we can extract components of xyz using the [] operator :
	FunctionOnMesh::Function & x = xyz[0], & y = xyz[1], & z = xyz[2];

	auto & sphere = R3.implicit ( x*x + y*y + z*z - 1. );
	double s2 = sqrt(0.5), s3 = sqrt(0.33333333333333);
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + *x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.5320889;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + 2*x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.5;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + 3*x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.47532976;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + 8*x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.401279;

	// Let's mesh an octant of a sphere
	// first, three corners, three middles of segments and the center :
	auto & A    = Cell::point();  x == 1.;   y == 0.;   z == 0.;
	auto & B    = Cell::point();  x == 0.;   y == 1.;   z == 0.;
	auto & C    = Cell::point();  x == 0.;   y == 0.;   z == 1.;
	auto & mAB  = Cell::point();  x == s2;   y == s2;   z == 0.;
	auto & mAC  = Cell::point();  x == s2;   y == 0.;   z == s2;
	auto & mBC  = Cell::point();  x == 0.;   y == s2;   z == s2;
	auto & mtri = Cell::point();  x == s3;   y == s3;   z == s3;

	// now build nine segments :
	auto & AmAB    = Mesh::segment ( A,   mAB,  10 );
	auto & AmAC    = Mesh::segment ( A,   mAC,  10 );
	auto & BmAB    = Mesh::segment ( B,   mAB,  10 );
	auto & BmBC    = Mesh::segment ( B,   mBC,  10 );
	auto & CmAC    = Mesh::segment ( C,   mAC,  10 );
	auto & CmBC    = Mesh::segment ( C,   mBC,  10 );
	auto & mABmtri = Mesh::segment ( mAB, mtri, 10 );
	auto & mACmtri = Mesh::segment ( mAC, mtri, 10 );
	auto & mBCmtri = Mesh::segment ( mBC, mtri, 10 );
	
	// now three rectangles :
	auto & rect_A = Mesh::rectangle ( mABmtri, mACmtri.reverse(), AmAC.reverse(), AmAB );
	auto & rect_B = Mesh::rectangle ( BmBC, mBCmtri, mABmtri.reverse(), BmAB.reverse() );
	auto & rect_C = Mesh::rectangle ( CmAC, mACmtri, mBCmtri.reverse(), CmBC.reverse() );

	// and finally join the rectangles :
	auto & octant_of_sphere = Mesh::join ( rect_A, rect_B, rect_C );
	
  octant_of_sphere.export_msh ("octant-of-sphere.msh");
	
	 cout << "reached end" << endl;
}
