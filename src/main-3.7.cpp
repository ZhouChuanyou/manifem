
// example presented in paragraph 3.7 of the manual
// meshes progressively a complicated surface, a kind of convolution between three tori

#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	cout << "this example takes some time" << endl;
	Function f1 = x*x + y*y + 0.1;
	Function f2 = 1. - power ( f1, -0.5 );
	Function d1 = z*z + f1 * f2 * f2;
	Function f3 = x*x + z*z + 0.1;
	Function f4 = 1. - power ( f3, -0.5 );
	Function d2 = y*y + f3 * f4 * f4;
	Function f5 = y*y + z*z + 0.1;
	Function f6 = 1. - power ( f5, -0.5 );
	Function d3 = x*x + f5 * f6 * f6;
	RR3.implicit ( smooth_min ( d1, d2, d3, tag::threshold, 0.1 ) == 0.1 );

	Mesh three_tori ( tag::progressive, tag::desired_length, 0.04 );

	three_tori.export_msh ("three-tori.msh");
	std::cout << "produced file three-tori.msh" << std::endl;
}
