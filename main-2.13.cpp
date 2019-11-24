
#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold spiral ( tag::Euclid, tag::of_dim, 1 );
	Function t = spiral.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	const double pi = 4.*atan(1.);
	
	Cell A ( tag::vertex );  t(A) = pi/2.;
	Cell B ( tag::vertex );  t(B) = 5.*pi;
	Mesh arc_of_spiral ( tag::segment, A, B, 50 );

	Function x = t*cos(t), y = t*sin(t);
	spiral.set_coordinates ( x && y );

	arc_of_spiral.draw_ps ("spiral.eps");
	arc_of_spiral.export_msh ("spiral.msh");

		cout << "reached end" << endl;
}
