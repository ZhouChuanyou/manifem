
// example presented in paragraph 3.16 of the manual
// union between a sphere and a cylinder

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;


int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	const double rs = 1.;    // radius of the sphere
	const double rc = 0.45;  // radius of the cylinder
	const double seg_size = 0.1;
	
	Manifold cylinder = RR3.implicit ( y*y + (z-0.5)*(z-0.5) == rc*rc );

	// base of the cylinder :
	cylinder.implicit ( x == 1.5 );

	Cell start1 ( tag::vertex );  x(start1) = 1.5;  y(start1) = 0.;  z(start1) = 0.5 + rc;
	std::vector < double > tau { 0., 1., 0. };
	Mesh circle_1 ( tag::progressive, tag::start_at, start1, tag::towards, tau,
	                tag::desired_length, seg_size                              );

	// intersection with the sphere
	Manifold intersection = cylinder.implicit ( x*x + y*y + z*z == rs*rs );
	Cell start2 ( tag::vertex );  x(start2) = 1.;  y(start2) = 0.;  z(start2) = 0.5 - rc;
	intersection.project ( start2 );
	tau = { 0., -1., 0. };
	Mesh circle_2 ( tag::progressive, tag::start_at, start2, tag::towards, tau,
	                tag::desired_length, seg_size                               );

	Mesh circles ( tag::join, circle_1, circle_2.reverse() );
	cylinder.set_as_working_manifold();
	tau = { -1., 0., 0. };
	Mesh cyl ( tag::progressive, tag::boundary, circles, tag::start_at, start1,
	           tag::towards, tau, tag::desired_length, seg_size                 );

	RR3.implicit ( x*x + y*y + z*z == rs*rs );
	tau = { 0., 0., -1. };
	Mesh sph ( tag::progressive, tag::boundary, circle_2,
	           tag::start_at, start2, tag::towards, tau,
	           tag::desired_length, seg_size              );

	Mesh all ( tag::join, cyl, sph );
	all.export_msh ("sphere-cyl.msh");
	cout << "produced file sphere-cyl.msh" << endl;
}
