
#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

  Manifold nut = RR3.implicit ( x*x + y*y + z*z + 1.5*x*y*z == 1. );

	// let's mesh a hemisphere (much deformed)
	// we take the four triangles in example 2.6 in the manual
	// and turn them into twelve rectangles
	// first, corners, middles of segments and centers of the four triangles :
	Cell S     ( tag::vertex );  x(S)     =  0.;  y(S)     = -1.;  z(S)     =  0.;
	Cell E     ( tag::vertex );  x(E)     =  1.;  y(E)     =  0.;  z(E)     =  0.;
	Cell N     ( tag::vertex );  x(N)     =  0.;  y(N)     =  1.;  z(N)     =  0.;
	Cell W     ( tag::vertex );  x(W)     = -1.;  y(W)     =  0.;  z(W)     =  0.;
	Cell up    ( tag::vertex );  x(up)    =  0.;  y(up)    =  0.;  z(up)    = 1.;
	// no need to project these
	Cell mSW   ( tag::vertex );  x(mSW)   = -1.;  y(mSW)   = -1.;  z(mSW)   =  0.;
	Cell mSE   ( tag::vertex );  x(mSE)   =  1.;  y(mSE)   = -1.;  z(mSE)   =  0.;
	Cell mNE   ( tag::vertex );  x(mNE)   =  1.;  y(mNE)   =  1.;  z(mNE)   =  0.;
	Cell mNW   ( tag::vertex );  x(mNW)   = -1.;  y(mNW)   =  1.;  z(mNW)   =  0.;
	nut.project ( mSW );  nut.project ( mSE );  nut.project ( mNE );  nut.project ( mNW );
	Cell mSup  ( tag::vertex );  x(mSup)  =  0.;  y(mSup)  = -1.;  z(mSup)  =  1.;
	Cell mEup  ( tag::vertex );  x(mEup)  =  1.;  y(mEup)  =  0.;  z(mEup)  =  1.;
	Cell mNup  ( tag::vertex );  x(mNup)  =  0.;  y(mNup)  =  1.;  z(mNup)  =  1.;
	Cell mWup  ( tag::vertex );  x(mWup)  = -1.;  y(mWup)  =  0.;  z(mWup)  =  1.;
	nut.project ( mSup );  nut.project ( mEup );  nut.project ( mNup );  nut.project ( mWup );
	Cell mSWup ( tag::vertex );  x(mSWup) = -1.;  y(mSWup) = -1.;  z(mSWup) =  1.;
	Cell mSEup ( tag::vertex );  x(mSEup) =  1.;  y(mSEup) = -1.;  z(mSEup) =  1.;
	Cell mNEup ( tag::vertex );  x(mNEup) =  1.;  y(mNEup) =  1.;  z(mNEup) =  1.;
	Cell mNWup ( tag::vertex );  x(mNWup) = -1.;  y(mNWup) =  1.;  z(mNWup) =  1.;
	nut.project ( mSWup );  nut.project ( mSEup );
	nut.project ( mNEup );  nut.project ( mNWup );

	// now build segments :
	int n = 10;
	Mesh S_mSW      ( tag::segment, S,    mSW,   n );
	Mesh S_mSE      ( tag::segment, S,    mSE,   n );
	Mesh E_mSE      ( tag::segment, E,    mSE,   n );
	Mesh E_mNE      ( tag::segment, E,    mNE,   n );
	Mesh N_mNE      ( tag::segment, N,    mNE,   n );
	Mesh N_mNW      ( tag::segment, N,    mNW,   n );
	Mesh W_mSW      ( tag::segment, W,    mSW,   n );
	Mesh W_mNW      ( tag::segment, W,    mNW,   n );
	Mesh S_mSup     ( tag::segment, S,    mSup,  n );
	Mesh N_mNup     ( tag::segment, N,    mNup,  n );
	Mesh E_mEup     ( tag::segment, E,    mEup,  n );
	Mesh W_mWup     ( tag::segment, W,    mWup,  n );
	Mesh mSW_mSWup  ( tag::segment, mSW,  mSWup, n );
	Mesh mSE_mSEup  ( tag::segment, mSE,  mSEup, n );
	Mesh mNE_mNEup  ( tag::segment, mNE,  mNEup, n );
	Mesh mNW_mNWup  ( tag::segment, mNW,  mNWup, n );
	Mesh up_mSup    ( tag::segment, up,   mSup,  n );
	Mesh up_mEup    ( tag::segment, up,   mEup,  n );
	Mesh up_mNup    ( tag::segment, up,   mNup,  n );
	Mesh up_mWup    ( tag::segment, up,   mWup,  n );
	Mesh mSup_mSEup ( tag::segment, mSup, mSEup, n );
	Mesh mSup_mSWup ( tag::segment, mSup, mSWup, n );
	Mesh mEup_mSEup ( tag::segment, mEup, mSEup, n );
	Mesh mEup_mNEup ( tag::segment, mEup, mNEup, n );
	Mesh mNup_mNEup ( tag::segment, mNup, mNEup, n );
	Mesh mNup_mNWup ( tag::segment, mNup, mNWup, n );
	Mesh mWup_mSWup ( tag::segment, mWup, mSWup, n );
	Mesh mWup_mNWup ( tag::segment, mWup, mNWup, n );
	
	// now the twelve rectangles :
	Mesh rect_S_SE  ( tag::rectangle,
                    S_mSE, mSE_mSEup, mSup_mSEup.reverse(), S_mSup.reverse() );
	Mesh rect_S_SW  ( tag::rectangle,
                    S_mSup, mSup_mSWup, mSW_mSWup.reverse(), S_mSW.reverse() );
	Mesh rect_E_SE  ( tag::rectangle,
                    E_mEup, mEup_mSEup, mSE_mSEup.reverse(), E_mSE.reverse() );
	Mesh rect_E_NE  ( tag::rectangle,
                    E_mNE, mNE_mNEup, mEup_mNEup.reverse(), E_mEup.reverse() );
	Mesh rect_N_NE  ( tag::rectangle,
                    N_mNup, mNup_mNEup, mNE_mNEup.reverse(), N_mNE.reverse() );
	Mesh rect_N_NW  ( tag::rectangle,
                    N_mNW, mNW_mNWup, mNup_mNWup.reverse(), N_mNup.reverse() );
	// yes, this sounds a lot like Carry Grant ...
	Mesh rect_W_SW  ( tag::rectangle,
                    W_mSW, mSW_mSWup, mWup_mSWup.reverse(), W_mWup.reverse() );
	Mesh rect_W_NW  ( tag::rectangle,
                     W_mWup, mWup_mNWup, mNW_mNWup.reverse(), W_mNW.reverse() );
	Mesh rect_up_SW ( tag::rectangle,
                    up_mWup, mWup_mSWup, mSup_mSWup.reverse(), up_mSup.reverse() );
	Mesh rect_up_SE ( tag::rectangle,
                    up_mSup, mSup_mSEup, mEup_mSEup.reverse(), up_mEup.reverse() );
	Mesh rect_up_NE ( tag::rectangle,
                    up_mEup, mEup_mNEup, mNup_mNEup.reverse(), up_mNup.reverse() );
	Mesh rect_up_NW ( tag::rectangle,
                    up_mNup, mNup_mNWup, mWup_mNWup.reverse(), up_mWup.reverse() );
											 
	// and finally join the rectangles :
	Mesh bumpy ( tag::join, list<Mesh>
		{ rect_S_SE, rect_S_SW, rect_E_SE, rect_E_NE, rect_N_NE, rect_N_NW,
		  rect_W_SW, rect_W_NW, rect_up_SW, rect_up_SE, rect_up_NE, rect_up_NW } );
	
  bumpy.export_msh ("bumpy.msh");
	
	 cout << "reached end" << endl;
}
