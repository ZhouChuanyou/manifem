
#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	Manifold sphere = RR3.implicit ( x*x + y*y + z*z == 1. );

	// let's mesh half of a sphere
	Cell E ( tag::vertex );  x(E) =  1.;   y(E) =  0.;   z(E) = 0.;
	Cell N ( tag::vertex );  x(N) =  0.;   y(N) =  1.;   z(N) = 0.;
	Cell W ( tag::vertex );  x(W) = -1.;   y(W) =  0.;   z(W) = 0.;
	Cell S ( tag::vertex );  x(S) =  0.;   y(S) = -1.;   z(S) = 0.;
	Cell up( tag::vertex );  x(up)=  0.;   y(up)=  0.;   z(up)= 1.;
	int n = 15;
	Mesh EN ( tag::segment, E, N, n );
	Mesh NW ( tag::segment, N, W, n );
	Mesh WS ( tag::segment, W, S, n );
	Mesh SE ( tag::segment, S, E, n );
	Mesh upE ( tag::segment, up, E, n );
	Mesh upN ( tag::segment, up, N, n );
	Mesh upW ( tag::segment, up, W, n );
	Mesh upS ( tag::segment, up, S, n );
	
	// now four triangles
	Mesh ENup ( tag::triangle, EN, upN.reverse(), upE );
	Mesh NWup ( tag::triangle, NW, upW.reverse(), upN );
	Mesh WSup ( tag::triangle, WS, upS.reverse(), upW );
	Mesh SEup ( tag::triangle, SE, upE.reverse(), upS );

	// and finally join the triangles :
	Mesh hemisphere ( tag::join, ENup, NWup, WSup, SEup );
	
  hemisphere.export_msh ("hemisphere.msh");
	
	cout << "reached end" << endl;
}
