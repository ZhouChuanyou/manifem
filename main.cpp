
#include "maniFEM.h"

using namespace maniFEM;

int main ()

{	// we choose our (topological) maximum dimension :
	Mesh::set_max_dim (2);
	// if we do not specify it, maniFEM will set topological dimension to 3

	// we choose our (geometric) space dimension :
	Manifold RR2 ( tag::Euclid, tag::of_dim, 2 );
	
	Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xy[0], y = xy[1];

	std::cout << (x/y==1).lhs->repr() << std::endl;
	
	Cell A ( tag::point );  A.core->name = "A";  x(A) = 1.; y(A) = -1.;
	Cell B ( tag::point );  B.core->name = "B";  x(B) = 1.; y(B) =  1.;
	Cell C ( tag::point );  C.core->name = "C";  x(C) = 0.; y(C) =  1.2;
	Cell D ( tag::point );  D.core->name = "D";  x(D) = 0.; y(D) = -1.;

	Mesh AB ( tag::segment, A, B, 10 );
	Mesh BC ( tag::segment, B, C, 15 );
	Mesh CD ( tag::segment, C, D, 10 );
	Mesh DA ( tag::segment, D, A, 15 );

	Mesh ABCD ( tag::pretty, tag::rectangle, AB, BC, CD, DA );

	ABCD.draw_ps ( "rectangle.eps");

	std::cout << "reached end" << std::endl;
}

