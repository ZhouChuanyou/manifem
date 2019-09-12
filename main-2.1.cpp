
#include "Mesh.h"

using namespace std;

int main () {

   auto RR2 = Manifold::euclid (2);  // Manifold RR2
   auto & xy = RR2.coordinate_system ("Lagrange degree one");
   // FunctionOnMesh::Function & xy
   auto & x = xy[0], & y = xy[1];  // FuncionOnMesh::Function & x, & y
   auto & A = Cell::point ("A"); x == -1; y == 0;   // Cell & A
   auto & C = Cell::point ("C"); x == 0; y == 0.5;
   auto & D = Cell::point ("D"); x == -1; y == 0.5;
   auto & E = Cell::point ("E"); x == 0; y == 1;
   auto & F = Cell::point ("F"); x == -1; y == 1;
   auto & G = Cell::point ("G"); x == 1; y == 0;
   auto & H = Cell::point ("H"); x == 1; y == 0.5;
   auto & AG = Mesh::segment ( A, G, 22 );   // Mesh & AG
   auto & GH = Mesh::segment ( G, H, 8 );
   auto & HC = Mesh::segment ( H, C, 12 );
   auto & CD = Mesh::segment ( C, D, 10 );
   auto & HD = Mesh::join ( HC, CD );        // Mesh & HD
   auto & DA = Mesh::segment ( D, A, 8 );
   auto & CE = Mesh::segment ( C, E, 7 );
   auto & EF = Mesh::segment ( E, F, 10 );
   auto & FD = Mesh::segment ( F, D, 7 );
   auto & AGHD = Mesh::rectangle ( AG, GH, HD, DA );   // Mesh & ABCD
   auto & CEFD = Mesh::rectangle ( CE, EF, FD, CD.reverse() );
   auto & L_shaped = Mesh::join ( AGHD, CEFD );   // Mesh & L_shaped
   // L_shaped.export_msh ("L-shaped.msh");
	 L_shaped.draw_ps ("L-shaped.eps");

	 cout << "reached end" << endl;
}
