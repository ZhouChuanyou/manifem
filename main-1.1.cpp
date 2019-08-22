
#include "Mesh.h"

using namespace std;

void do_test ();
// Cell * Cell::celula_especial;

int main () {

	// we choose our (geometric) space dimension :
	auto & environment = Manifold::euclid (3);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'environment' :
	auto & xyz = environment.coordinate_system ("Lagrange degree one");

	// we can extract components of xyz using the [] operator :
	auto & x = xyz[0], & y = xyz[1], & z = xyz[2];

	x.components[0]->name = "x";
	y.components[0]->name = "y";
	z.components[0]->name = "z";

	// Let's build a rectangular mesh. First, the four corners :
	auto & SW = Cell::point ("SW");  x ==  0.;  y ==  0.;  z ==  0.;
	auto & SE = Cell::point ("SE");  x ==  1.;  y ==  0.;  z ==  0.;
	auto & NE = Cell::point ("NE");  x ==  1.;  y ==  1.;  z ==  1.;
	auto & NW = Cell::point ("NW");  x ==  0.;  y ==  1.;  z ==  0.;

	// we access the coordinates of a point using the () operator :
	cout << "coordinates of fourth corner : " << x(NW) << " " << y(NW) << " " << z(NW) << endl;
	
	// Now build the four sides of the rectangle :
	auto & south = Mesh::segment ( SW, SE, 10 );
	auto & east  = Mesh::segment ( SE, NE, 10 );
	auto & north = Mesh::segment ( NE, NW, 10 );
	auto & west  = Mesh::segment ( NW, SW, 10 );
	// And now the rectangle :
	auto & rect_mesh = Mesh::rectangle ( south, east, north, west );
	 
	// We may want to visualize the resulting mesh.
	// Here is one way to export the mesh in the "msh" format :
	rect_mesh.export_msh ("rectangle.msh");
	// rect_mesh.draw_ps ("rectangle.eps");
	
	// Let's define a symbolic function to integrate
	auto & f = x*x+1/(5+y), & g = x*y;

	// and compute its integral on the rectangle, using Gauss quadrature with 9 points :
	auto & integ = Integrator::gauss ("Q9");
	//	Cell::celula_especial = & NW;
	cout << "integral of " << f.repr() << " = " << f.integrate ( rect_mesh, integ ) << endl;
	cout << "integral of " << g.repr() << " = " << g.integrate ( rect_mesh, integ ) << endl;
	
	// new points, segments and rectangles are created by these calls
	// the fields 'x', 'y' and 'z' are interpolated in the new vertices
	// by using the information stored in the function xyz

	// abstract variational formulation :
	// auto & uu = FunctionOnMesh::unknown ( u, "Lagrange degree one");
	// The values of 'uu' are related to the field 'u'.
	// On segments, 'uu' varies linearly. On a triangle, it varies
	// also linearly. On a rectangle, it is a polynomial of degree 1,
	// so it has a linear part and a bilinear one.
	// auto & w = FunctionOnMesh::test ( uu );
	// 'w' will be used for defining variational formulations in tandem with 'uu'

	// auto & integ = Integrator::gauss ("Q9");
	// integ nao recebe as coordenadas espaciais x e y
	// fe (mais abaixo) extrai essa informacao da formulacao variacional
	// cout << ( xx*yy ) .integrate ( malha, integ ) << endl;
	
/*
	auto & var_pb =
		( uu.deriv(xx)*w.deriv(xx) + uu.deriv(yy)*w.deriv(yy) ) .integrate(malha)
		== w.integrate(malha) + w.integrate(south);

	var_pb.prescribe_on (north);  uu == 0.;          w == 0.;
	var_pb.prescribe_on (east);   uu == xx*(1.-yy);  w == 0.;
	var_pb.prescribe_on (west);   uu == 0.;          w == 0.;

	auto & fe = FiniteElement::lagrange ("Q1", malha );
	// fe.set_integrator ("on cells of dimension", 1, "gauss", 3, "nodes");
	fe.set_integrator ("on cells of dimension", 2, "gauss", 9, "nodes");
	var_pb.set_finite_element ( fe );
	// A finite element should perform the following tasks.
	// Build a base in the (discrete) Hilbert space and index it;
	// the index should be somehow related to an enumeration of cells.
	// Take a variational problem and express it as a system of linear equations,
	// by replacing the unknown and the test functions by functions in the base
	// and by computing the respective integrals.
	// Take the solution of the system and interpret it as a function on the mesh.

	var_pb.discretize ();
*/
//	do_test ();
	
	 cout << "reached end" << endl;
}
