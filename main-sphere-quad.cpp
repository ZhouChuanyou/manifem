
#include "Mesh.h"
#include "math.h"

using namespace std;

int main () {

	// we choose our (geometric) space dimension :
	auto RR3 = Manifold::euclid (3);
	// for now, the constructor Manifold() sets Mesh::environment to 'this'
	// this may change in the future
	
	// xyz is a map defined on our future mesh with values in 'R3' :
	auto & xyz = RR3.coordinate_system ("Lagrange degree one");
	xyz.components[0]->name = "x";
	xyz.components[1]->name = "y";
	xyz.components[2]->name = "z";

	// we can extract components of xyz using the [] operator :
	FunctionOnMesh::Function & x = xyz[0], & y = xyz[1], & z = xyz[2];

	// auto & sphere = R3.implicit ( x*x + y*y + z*z - 1. );
	// double s2 = sqrt(0.5), s3 = sqrt(0.33333333333333), s4 = s3;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.5320889, s4 = 0.652704;
	auto & sphere = R3.implicit ( x*x + y*y + z*z + 1.5*x*y*z - 1. );
	double s2 = sqrt(0.5), s3 = 0.5148688, s4 = 0.7223518;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + 1.8*x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.505704, s4 = 0.8011915;
	// auto & sphere = R3.implicit ( x*x + y*y + z*z + 2.*x*y*z - 1. );
	// double s2 = sqrt(0.5), s3 = 0.5, s4 = 1.;

	// let's mesh a hemisphere (much deformed)
	// first, three corners, three middles of segments and the center :
	auto & A    = Cell::point();  x ==  1.;   y ==  0.;   z == 0.;
	auto & B    = Cell::point();  x ==  0.;   y ==  1.;   z == 0.;
	auto & C    = Cell::point();  x == -1.;   y ==  0.;   z == 0.;
	auto & D    = Cell::point();  x ==  0.;   y == -1.;   z == 0.;
	auto & N    = Cell::point();  x ==  0.;   y ==  0.;   z == 1.;
	auto & mAB  = Cell::point();  x ==  s2;   y ==  s2;   z == 0.;
	auto & mAN  = Cell::point();  x ==  s2;   y ==  0.;   z == s2;
	auto & mBN  = Cell::point();  x ==  0.;   y ==  s2;   z == s2;
	auto & mABN = Cell::point();  x ==  s3;   y ==  s3;   z == s3;
	auto & mBC  = Cell::point();  x == -s2;   y ==  s2;   z == 0.;
	auto & mCN  = Cell::point();  x == -s2;   y ==  0.;   z == s2;
	auto & mBCN = Cell::point();  x == -s4;   y ==  s4;   z == s4;
	auto & mCD  = Cell::point();  x == -s2;   y == -s2;   z == 0.;
	auto & mDN  = Cell::point();  x ==  0.;   y == -s2;   z == s2;
	auto & mCDN = Cell::point();  x == -s3;   y == -s3;   z == s3;
	auto & mDA  = Cell::point();  x ==  s2;   y == -s2;   z == 0.;
	auto & mDAN = Cell::point();  x ==  s4;   y == -s4;   z == s4;

	// now build segments :
	int n = 10;
	auto & AmAB    = Mesh::segment ( A,   mAB,  n );
	auto & AmAN    = Mesh::segment ( A,   mAN,  n );
	auto & BmAB    = Mesh::segment ( B,   mAB,  n );
	auto & BmBN    = Mesh::segment ( B,   mBN,  n );
	auto & NmAN    = Mesh::segment ( N,   mAN,  n );
	auto & NmBN    = Mesh::segment ( N,   mBN,  n );
	auto & mABmABN = Mesh::segment ( mAB, mABN, n );
	auto & mANmABN = Mesh::segment ( mAN, mABN, n );
	auto & mBNmABN = Mesh::segment ( mBN, mABN, n );
	auto & BmBC    = Mesh::segment ( B,   mBC,  n );
	auto & CmBC    = Mesh::segment ( C,   mBC,  n );
	auto & CmCN    = Mesh::segment ( C,   mCN,  n );
	auto & NmCN    = Mesh::segment ( N,   mCN,  n );
	auto & mBCmBCN = Mesh::segment ( mBC, mBCN, n );
	auto & mBNmBCN = Mesh::segment ( mBN, mBCN, n );
	auto & mCNmBCN = Mesh::segment ( mCN, mBCN, n );
	auto & CmCD    = Mesh::segment ( C,   mCD,  n );
	auto & DmCD    = Mesh::segment ( D,   mCD,  n );
	auto & DmDN    = Mesh::segment ( D,   mDN,  n );
	auto & NmDN    = Mesh::segment ( N,   mDN,  n );
	auto & mCDmCDN = Mesh::segment ( mCD, mCDN, n );
	auto & mCNmCDN = Mesh::segment ( mCN, mCDN, n );
	auto & mDNmCDN = Mesh::segment ( mDN, mCDN, n );
	auto & DmDA    = Mesh::segment ( D,   mDA,  n );
	auto & AmDA    = Mesh::segment ( A,   mDA,  n );
	auto & mDAmDAN = Mesh::segment ( mDA, mDAN, n );
	auto & mDNmDAN = Mesh::segment ( mDN, mDAN, n );
	auto & mANmDAN = Mesh::segment ( mAN, mDAN, n );
	
	// now three rectangles in each octant :
	auto & rect_A_to_B  = Mesh::rectangle ( mABmABN, mANmABN.reverse(), AmAN.reverse(), AmAB );
	auto & rect_B_to_A  = Mesh::rectangle ( mBNmABN, mABmABN.reverse(), BmAB.reverse(), BmBN );
	auto & rect_N_to_AB = Mesh::rectangle ( mANmABN, mBNmABN.reverse(), NmBN.reverse(), NmAN );
	auto & rect_B_to_C  = Mesh::rectangle ( mBCmBCN, mBNmBCN.reverse(), BmBN.reverse(), BmBC );
	auto & rect_C_to_B  = Mesh::rectangle ( mCNmBCN, mBCmBCN.reverse(), CmBC.reverse(), CmCN );
	auto & rect_N_to_BC = Mesh::rectangle ( mBNmBCN, mCNmBCN.reverse(), NmCN.reverse(), NmBN );
	auto & rect_C_to_D  = Mesh::rectangle ( mCDmCDN, mCNmCDN.reverse(), CmCN.reverse(), CmCD );
	auto & rect_D_to_C  = Mesh::rectangle ( mDNmCDN, mCDmCDN.reverse(), DmCD.reverse(), DmDN );
	auto & rect_N_to_CD = Mesh::rectangle ( mCNmCDN, mDNmCDN.reverse(), NmDN.reverse(), NmCN );
	auto & rect_D_to_A  = Mesh::rectangle ( mDAmDAN, mDNmDAN.reverse(), DmDN.reverse(), DmDA );
	auto & rect_A_to_D  = Mesh::rectangle ( mANmDAN, mDAmDAN.reverse(), AmDA.reverse(), AmAN );
	auto & rect_N_to_DA = Mesh::rectangle ( mDNmDAN, mANmDAN.reverse(), NmAN.reverse(), NmDN );

	// and finally join the rectangles :
	auto & hemisphere = Mesh::join ( list<Mesh*>{ &rect_A_to_B, &rect_B_to_A, &rect_N_to_AB, &rect_B_to_C, &rect_C_to_B, &rect_N_to_BC, &rect_C_to_D, &rect_D_to_C, &rect_N_to_CD, &rect_D_to_A, &rect_A_to_D, &rect_N_to_DA } );
	
  hemisphere.export_msh ("hemisphere.msh");
	
	 cout << "reached end" << endl;
}
