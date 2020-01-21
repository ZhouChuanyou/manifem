
// example presented in paragraph 2.2 of the manual
// builds an L-shaped mesh which mixed traingles and rectangles

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
	Mesh AB ( tag::segment, A, B, tag::divided_in, 10 );
	Mesh BC ( tag::segment, B, C, tag::divided_in,  8 );
	Mesh CD ( tag::segment, C, D, tag::divided_in, 10 );
	Mesh DA ( tag::segment, D, A, tag::divided_in,  8 );
	Mesh CE ( tag::segment, C, E, tag::divided_in,  7 );
	Mesh EF ( tag::segment, E, F, tag::divided_in, 10 );
	Mesh FD ( tag::segment, F, D, tag::divided_in,  7 );
	Mesh BG ( tag::segment, B, G, tag::divided_in, 12 );
	Mesh GH ( tag::segment, G, H, tag::divided_in,  8 );
	Mesh HC ( tag::segment, H, C, tag::divided_in, 12 );

	Mesh ABCD ( tag::rectangle, AB, BC, CD, DA );
	Mesh CEFD ( tag::rectangle, CE, EF, FD, CD.reverse() );
	Mesh BGHC ( tag::rectangle, BG, GH, HC, BC.reverse(), tag::with_triangles );
	Mesh L_shaped ( tag::join, ABCD, CEFD, BGHC );

	L_shaped.draw_ps ( "L-shaped.eps");
	L_shaped.export_msh ("L-shaped.msh");
	
	cout << "produced files L-shaped.eps and L-shaped.msh" << endl;
}