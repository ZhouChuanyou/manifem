
#include "Mesh.h"

using namespace std;

int main () {

   auto RR3 = Manifold::euclid (2);  // Manifold RR3
   auto & xy = RR3.coordinate_system ("Lagrange degree one");
   // FunctionOnMesh::Function & xy
   auto & x = xy[0], & y = xy[1];  // FuncionOnMesh::Function & x, & y
   auto & A = Cell::point ("A");  x == -1. ;  y == 0.;   // Cell & A
   auto & B = Cell::point ("B");  x ==  0. ;  y == 0.;
   auto & C = Cell::point ("C");  x ==  1. ;  y == 0.;
   auto & D = Cell::point ("D");  x == -0.5;  y == 1.;
   auto & E = Cell::point ("E");  x ==  0.5;  y == 1.;
   auto & AB = Mesh::segment ( A, B, 8 );   // Mesh & AB
   auto & BC = Mesh::segment ( B, C, 8 );
   auto & AD = Mesh::segment ( A, D, 8 );
   auto & BD = Mesh::segment ( B, D, 8 );
   auto & BE = Mesh::segment ( B, E, 8 );
   auto & CE = Mesh::segment ( C, E, 8 );
   auto & ED = Mesh::segment ( E, D, 8 );
   auto & ABD = Mesh::triangle ( AB, BD, AD.reverse() );   // Mesh & ABD
   auto & BCE = Mesh::triangle ( BC, CE, BE.reverse() );
   auto & BED = Mesh::triangle ( BE, ED, BD.reverse() );
   auto & three_tri = Mesh::join ( ABD, BCE, BED );

   three_tri.export_msh ("three-tri.msh");
	
   cout << "reached end" << endl;
}
