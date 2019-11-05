
#include "maniFEM.h"
#include "math.h"

using namespace maniFEM;
using namespace std;

int main ()

{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

  Manifold nut = RR3.implicit ( x*x + y*y + z*z + 1.5*x*y*z == 1. );
	double s2 = sqrt(0.5), s3 = 0.5148688, s4 = 0.7223518;
	// auto & nut = RR3.implicit ( x*x + y*y + z*z + 1.8*x*y*z == 1. );
	// double s2 = sqrt(0.5), s3 = 0.505704, s4 = 0.8011915;
	// auto & nut = RR3.implicit ( x*x + y*y + z*z + 2.0*x*y*z == 1. );
	// double s2 = sqrt(0.5), s3 = 0.5, s4 = 1.;

	// let's mesh a hemisphere (much deformed)
	// we take the four triangles in example 2.6 in the manual
	// and split each of them in four quadranges
	// first, corners, middles of segments and centers of the four triangles :
	Cell A    ( tag::vertex );  x(A) =  1.;     y(A) =  0.;     z(A) = 0.;
	Cell B    ( tag::vertex );  x(B) =  0.;     y(B) =  1.;     z(B) = 0.;
	Cell C    ( tag::vertex );  x(C) = -1.;     y(C) =  0.;     z(C) = 0.;
	Cell D    ( tag::vertex );  x(D) =  0.;     y(D) = -1.;     z(D) = 0.;
	Cell N    ( tag::vertex );  x(N) =  0.;     y(D) =  0.;     z(D) = 1.;
	Cell mAB  ( tag::vertex );  x(mAB) =   s2;  y(mAB) =   s2;  z(mAB) =  0.;
	Cell mAN  ( tag::vertex );  x(mAN) =   s2;  y(mAN) =   0.;  z(mAN) =  s2;
	Cell mBN  ( tag::vertex );  x(mBN) =   0.;  y(mBN) =   s2;  z(mBN) =  s2;
	Cell mABN ( tag::vertex );  x(mABN) =  s3;  y(mABN) =  s3;  z(mABN) = s3;
	Cell mBC  ( tag::vertex );  x(mBC) =  -s2;  y(mBC) =   s2;  z(mBC) =  0.;
	Cell mCN  ( tag::vertex );  x(mCN) =  -s2;  y(mCN) =   0.;  z(mCN) =  s2;
	Cell mBCN ( tag::vertex );  x(mBCN) = -s4;  y(mBCN) =  s4;  z(mBCN) = s4;
	Cell mCD  ( tag::vertex );  x(mCD) =  -s2;  y(mCD) =  -s2;  z(mCD) =  0.;
	Cell mDN  ( tag::vertex );  x(mDN) =   0.;  y(mDN) =  -s2;  z(mDN) =  s2;
	Cell mCDN ( tag::vertex );  x(mCDN) = -s3;  y(mCDN) = -s3;  z(mCDN) = s3;
	Cell mDA  ( tag::vertex );  x(mDA) =   s2;  y(mDA) =  -s2;  z(mDA) =  0.;
	Cell mDAN ( tag::vertex );  x(mDAN) =  s4;  y(mDAN) = -s4;  z(mDAN) = s4;

	// now build segments :
	int n = 10;
	Mesh AmAB    ( tag::segment, A,   mAB,  n );
	Mesh AmAN    ( tag::segment, A,   mAN,  n );
	Mesh BmAB    ( tag::segment, B,   mAB,  n );
	Mesh BmBN    ( tag::segment, B,   mBN,  n );
	Mesh NmAN    ( tag::segment, N,   mAN,  n );
	Mesh NmBN    ( tag::segment, N,   mBN,  n );
	Mesh mABmABN ( tag::segment, mAB, mABN, n );
	Mesh mANmABN ( tag::segment, mAN, mABN, n );
	Mesh mBNmABN ( tag::segment, mBN, mABN, n );
	Mesh BmBC    ( tag::segment, B,   mBC,  n );
	Mesh CmBC    ( tag::segment, C,   mBC,  n );
	Mesh CmCN    ( tag::segment, C,   mCN,  n );
	Mesh NmCN    ( tag::segment, N,   mCN,  n );
	Mesh mBCmBCN ( tag::segment, mBC, mBCN, n );
	Mesh mBNmBCN ( tag::segment, mBN, mBCN, n );
	Mesh mCNmBCN ( tag::segment, mCN, mBCN, n );
	Mesh CmCD    ( tag::segment, C,   mCD,  n );
	Mesh DmCD    ( tag::segment, D,   mCD,  n );
	Mesh DmDN    ( tag::segment, D,   mDN,  n );
	Mesh NmDN    ( tag::segment, N,   mDN,  n );
	Mesh mCDmCDN ( tag::segment, mCD, mCDN, n );
	Mesh mCNmCDN ( tag::segment, mCN, mCDN, n );
	Mesh mDNmCDN ( tag::segment, mDN, mCDN, n );
	Mesh DmDA    ( tag::segment, D,   mDA,  n );
	Mesh AmDA    ( tag::segment, A,   mDA,  n );
	Mesh mDAmDAN ( tag::segment, mDA, mDAN, n );
	Mesh mDNmDAN ( tag::segment, mDN, mDAN, n );
	Mesh mANmDAN ( tag::segment, mAN, mDAN, n );
	
	// now three rectangles in each octant :
	Mesh rect_A_to_B  ( tag::rectangle, mABmABN, mANmABN.reverse(), AmAN.reverse(), AmAB );
	Mesh rect_B_to_A  ( tag::rectangle, mBNmABN, mABmABN.reverse(), BmAB.reverse(), BmBN );
	Mesh rect_N_to_AB ( tag::rectangle, mANmABN, mBNmABN.reverse(), NmBN.reverse(), NmAN );
	Mesh rect_B_to_C  ( tag::rectangle, mBCmBCN, mBNmBCN.reverse(), BmBN.reverse(), BmBC );
	Mesh rect_C_to_B  ( tag::rectangle, mCNmBCN, mBCmBCN.reverse(), CmBC.reverse(), CmCN );
	Mesh rect_N_to_BC ( tag::rectangle, mBNmBCN, mCNmBCN.reverse(), NmCN.reverse(), NmBN );
	Mesh rect_C_to_D  ( tag::rectangle, mCDmCDN, mCNmCDN.reverse(), CmCN.reverse(), CmCD );
	Mesh rect_D_to_C  ( tag::rectangle, mDNmCDN, mCDmCDN.reverse(), DmCD.reverse(), DmDN );
	Mesh rect_N_to_CD ( tag::rectangle, mCNmCDN, mDNmCDN.reverse(), NmDN.reverse(), NmCN );
	Mesh rect_D_to_A  ( tag::rectangle, mDAmDAN, mDNmDAN.reverse(), DmDN.reverse(), DmDA );
	Mesh rect_A_to_D  ( tag::rectangle, mANmDAN, mDAmDAN.reverse(), AmDA.reverse(), AmAN );
	Mesh rect_N_to_DA ( tag::rectangle, mDNmDAN, mANmDAN.reverse(), NmAN.reverse(), NmDN );

	// and finally join the rectangles :
	Mesh bumpy ( tag::join, list<Mesh>
	  { rect_A_to_B, rect_B_to_A, rect_N_to_AB, rect_B_to_C, rect_C_to_B,
	    rect_N_to_BC, rect_C_to_D, rect_D_to_C, rect_N_to_CD, rect_D_to_A,
	    rect_A_to_D, rect_N_to_DA                                           } );
	
  bumpy.export_msh ("bumpy.msh");
	
	 cout << "reached end" << endl;
}
