
#include "maniFEM.h"

using namespace maniFEM;
using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	
	// xy is a map defined on our future mesh with values in RR2 :
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

	// we can extract components of xyz using the [] operator :
	Function x = xy[0],  y = xy[1];

	Cell A ( tag::vertex );  x(A) = -1.;  y(A) = 0.;
	Cell B ( tag::vertex );  x(B) =  0.;  y(B) = 0.;
	Cell C ( tag::vertex );  x(C) =  0.;  y(C) = 0.5;
	Cell D ( tag::vertex );  x(D) = -1.;  y(D) = 0.5;
	Cell E ( tag::vertex );  x(E) =  0.;  y(E) = 1.;
	Cell F ( tag::vertex );  x(F) = -1.;  y(F) = 1.;
	Cell G ( tag::vertex );  x(G) =  1.;  y(G) = 0.;
	Cell H ( tag::vertex );  x(H) =  1.;  y(H) = 0.5;
	Mesh AB ( tag::segment, A, B, 10 );
  Mesh BC ( tag::segment, B, C, 8 );
	Mesh CD ( tag::segment, C, D, 10 );
  Mesh DA ( tag::segment, D, A, 8 );
	Mesh CE ( tag::segment, C, E, 7 );
	Mesh EF ( tag::segment, E, F, 10 );
	Mesh FD ( tag::segment, F, D, 7 );
	Mesh BG ( tag::segment, B, G, 12 );
	Mesh GH ( tag::segment, G, H, 8 );
	Mesh HC ( tag::segment, H, C, 12 );

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh CEFD ( tag::rectangle, CE, EF, FD, CD.reverse() );
	Mesh BGHC ( tag::rectangle, GH, HC, BC.reverse(), BG );
	Mesh L_shaped ( tag::join, ABCD, CEFD, BGHC );

  L_shaped.export_msh ("L-shaped.msh");
	L_shaped.draw_ps ( "L-shaped.eps");
	
	cout << "reached end" << endl;
}
