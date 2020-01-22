
// example presented in paragraph 3.22 of the manual
// the code shown in the manual does not work (yet)
// we fake the desired result by building a mesh on a paraboloid

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function  x = xyz[0],  y = xyz[1],  z = xyz[2];

	Manifold parab_manif = RR3.implicit ( z == std::sqrt(0.75) * power ( x*x + y*y, 1.4 ) );

	parab_manif.implicit ( x*x + y*y == 1. );
	Cell A ( tag::vertex );  x(A) = 1.;  y(A) = 0.;  z(A) = std::sqrt(0.75);
	Mesh circle ( tag::progressive, tag::start_at, A, tag::desired_length, 0.2 );

	parab_manif.set_as_working_manifold();
	std::vector < double > tau { 0., 0., -1. };
	Mesh parab ( tag::progressive, tag::boundary, circle,
							 tag::start_at, A, tag::towards, tau,
	             tag::desired_length, 0.17                 );

	RR3.set_coordinates ( x && y );
	parab.export_msh ("parab.msh");
	cout << "produced file parab.msh" << endl;
}
