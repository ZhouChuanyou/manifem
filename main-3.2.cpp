
#include "Mesh.h"
#include <cmath>

using namespace std;

// builds a ring-shaped mesh

int main () {

	Manifold RR2 = Manifold::euclid (2);
	FunctionOnMesh::Function & xy = RR2.coordinate_system ("Lagrange degree one");
	FunctionOnMesh::Function & x = xy[0], & y = xy[1];

	short int n_sectors = 15;
	double step_theta = 8*atan(1.)/n_sectors;
	short int radial_divisions = 10;
	short int rot_divisions = 5;
	
	// start the process by building a segment
	Cell & A = Cell::point (); x == 1.; y == 0.;
	Cell & B = Cell::point (); x == 2.; y == 0.;
	Mesh & ini_seg = Mesh::segment ( A, B, radial_divisions );
	// we will need to vary ini_seg, A and B, so we use pointers :
	Mesh * prev_seg_p = & ini_seg;
	Cell * A_p = &A, * B_p = &B;
	// 'l' will be a list of trapezoidal meshes
	list <Mesh*> l;

	for ( short int i = 1; i < n_sectors; i++ )
	{	double theta = i * step_theta;
		// we build two new points
		Cell & C = Cell::point (); x == cos(theta);    y == sin(theta);
		Cell & D = Cell::point (); x == 2.*cos(theta); y == 2.*sin(theta);
		// and three new segments
		Mesh & BD = Mesh::segment ( *B_p, D, rot_divisions );
		Mesh & DC = Mesh::segment ( D, C, radial_divisions );
		Mesh & CA = Mesh::segment ( C, *A_p, rot_divisions );
		Mesh & rect = Mesh::rectangle ( *prev_seg_p, BD, DC, CA );
		l.push_back ( &rect );
		prev_seg_p = & ( DC.reverse() );
		A_p = &C; B_p = &D;                                                  	}

	// let's build the last sector, and close the ring
	// prev_seg_p, A_p and B_p have rotated during the construction process
	// but ini_seg, A and B are the same, initial, ones
	Mesh & outer = Mesh::segment ( *B_p, B, rot_divisions );
	Mesh & inner = Mesh::segment ( A, *A_p, rot_divisions );
	Mesh & rect = Mesh::rectangle ( outer, ini_seg.reverse(), inner, *prev_seg_p );
	l.push_back ( &rect );
	
	Mesh & ring = Mesh::join ( l );
	ring.export_msh ("ring.msh");

	cout << "reached end" << endl;
}
