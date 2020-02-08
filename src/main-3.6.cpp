
// example presented in paragraph 3.6 of the manual
// meshes a sphere progressively

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function  x = xyz[0],  y = xyz[1],  z = xyz[2];

	cout << "this example takes some time" << endl;

	RR3.implicit ( x*x + y*y + z*z == 0.01 );
	Mesh sphere ( tag::progressive, tag::desired_length, 0.01 );

	sphere.export_msh ("sphere.msh");

	cout << "produced file sphere.msh" << endl;
}
