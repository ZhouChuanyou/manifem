
#include "Mesh.h"

using namespace std;

int main () {

   auto RR3 = Manifold::euclid (2);  // Manifold RR3
   auto & xy = RR3.coordinate_system ("Lagrange degree one");
   // FunctionOnMesh::Function & xy
   auto & x = xy[0], & y = xy[1];  // FuncionOnMesh::Function & x, & y
   auto & A = Cell::point ("A");  x == -1. ;  y == 0. ;   // Cell & A
   auto & B = Cell::point ("B");  x ==  0. ;  y == 0. ;
   auto & C = Cell::point ("C");  x ==  1. ;  y == 0. ;
   auto & D = Cell::point ("D");  x == -0.5;  y == 0.8;
   auto & E = Cell::point ("E");  x ==  0.5;  y == 0.8;
   auto & F = Cell::point ("F");  x ==  0. ;  y == 1. ;
   auto & AB = Mesh::segment ( A, B, 8 );   // Mesh & AB
   auto & BC = Mesh::segment ( B, C, 8 );
   auto & AD = Mesh::segment ( A, D, 8 );
   auto & BD = Mesh::segment ( B, D, 8 );
   auto & BE = Mesh::segment ( B, E, 8 );
   auto & CE = Mesh::segment ( C, E, 8 );
   auto & EF = Mesh::segment ( E, F, 8 );
   auto & FD = Mesh::segment ( F, D, 8 );
   auto & ABD = Mesh::triangle ( AB, BD, AD.reverse() );   // Mesh & ABD
   auto & BCE = Mesh::triangle ( BC, CE, BE.reverse() );
   auto & BEFD = Mesh::rectangle ( BE, EF, FD, BD.reverse() );
   auto & two_tri_one_rect = Mesh::join ( ABD, BEFD, BCE );

   two_tri_one_rect.export_msh ("two-tri-one-rect.msh");
	
   cout << "reached end" << endl;
}
