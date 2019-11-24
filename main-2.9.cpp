
#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0],  y = xy[1];

	Manifold manif_circle = RR2.implicit ( x*x + y*y == 1. );

	Cell N ( tag::vertex );  x(N) =  0.;  y(N) =  1.;
	Cell W ( tag::vertex );  x(W) = -1.;  y(W) =  0.;
  Cell S ( tag::vertex );  x(S) =  0.;  y(S) = -1.;
  Cell E ( tag::vertex );  x(E) =  1.;  y(E) =  0.;

	Mesh NW ( tag::segment, N, W, 10 );
	Mesh WS ( tag::segment, W, S, 10 );
	Mesh SE ( tag::segment, S, E, 10 );
	Mesh EN ( tag::segment, E, N, 10 );

	RR2.set_as_working_manifold();

	Mesh disk ( tag::rectangle, NW, WS, SE, EN );

	disk.draw_ps ("disk.eps");
  disk.export_msh ("disk.msh");
	
	 cout << "reached end" << endl;
}
