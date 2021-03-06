
// finite-elem.cpp 2021.02.27

//   This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//   Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com
//   https://github.com/cristian-barbarosie/manifem

//   ManiFEM is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, either version 3 of the License
//   or (at your option) any later version.

//   ManiFEM is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Lesser General Public License for more details.

//   You should have received a copy of the GNU Lesser General Public License
//   along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

#include "math.h"
#include <sstream>

#include "iterator.h"
#include "function.h"
#include "finite-elem.h"

using namespace maniFEM;

//-----------------------------------------------------------------------------------------//

Integrator::Gauss::Gauss ( const tag::gauss_quadrature & q,
                           const tag::FromFiniteElementWithMaster &, FiniteElement & fe )

: Integrator::Core ()

// http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf  ../fea6.pdf
// E.B. Becker, G.F. Carey, J.T. Oden, Finite Elements, an introduction, vol 1

{ FiniteElement::WithMaster * fe_core = dynamic_cast < FiniteElement::WithMaster * > ( fe.core );
	assert ( fe_core );
	Manifold master_manifold = fe_core->master_manif;
	Function xi_eta = master_manifold.coordinates();
	assert ( xi_eta.nb_of_components() == 2 );
	Function xi = xi_eta[0], eta = xi_eta[1];
	switch ( q )
	{	case tag::tri_3 :  // Gauss quadrature with three points on a triangle
		{	FiniteElement::WithMaster::Triangle * tri_fe =
				dynamic_cast < FiniteElement::WithMaster::Triangle * > ( fe.core );
			assert ( tri_fe );
			this->points.reserve(3);  this->weights.reserve(3);
			// points
			Cell Gauss_3_O ( tag::vertex );  xi(Gauss_3_O) = 1./6.;  eta(Gauss_3_O) = 1./6.;
			Cell Gauss_3_A ( tag::vertex );  xi(Gauss_3_A) = 2./3.;  eta(Gauss_3_A) = 1./6.;
			Cell Gauss_3_B ( tag::vertex );  xi(Gauss_3_B) = 1./6.;  eta(Gauss_3_B) = 2./3.;
			// weights
			double Gw_3_O = 1./6., Gw_3_A = 1./6., Gw_3_B = 1./6.;
			// in the book of Oden we have the same weights but a different distribution of points
			// instead of  2/3 1/6 1/6  Oden has  middle of segments 0 0.5 0.5
			// maybe the points are not unique ?
			this->points.push_back ( Gauss_3_O );
			this->weights.push_back ( Gw_3_O );
			this->points.push_back ( Gauss_3_A );
			this->weights.push_back ( Gw_3_A );
			this->points.push_back ( Gauss_3_B );
			this->weights.push_back ( Gw_3_B );
			break;                                                                             }
		case tag::tri_3_Oden :
		// Gauss quadrature with three points on a triangle, points presented in the book
		// E.B. Becker, G.F. Carey, J.T. Oden, Finite Elements, an introduction, vol 1
		{	FiniteElement::WithMaster::Triangle * tri_fe =
				dynamic_cast < FiniteElement::WithMaster::Triangle * > ( fe.core );
			assert ( tri_fe );
			this->points.reserve(3);  this->weights.reserve(3);
			// points
			Cell Gauss_3_O ( tag::vertex );  xi(Gauss_3_O) = 0.5;  eta(Gauss_3_O) = 0.5;
			Cell Gauss_3_A ( tag::vertex );  xi(Gauss_3_A) = 0.;   eta(Gauss_3_A) = 0.5;
			Cell Gauss_3_B ( tag::vertex );  xi(Gauss_3_B) = 0.5;  eta(Gauss_3_B) = 0.;
			// weights
			double Gw_3_O = 1./6., Gw_3_A = 1./6., Gw_3_B = 1./6.;
			this->points.push_back ( Gauss_3_O );
			this->weights.push_back ( Gw_3_O );
			this->points.push_back ( Gauss_3_A );
			this->weights.push_back ( Gw_3_A );
			this->points.push_back ( Gauss_3_B );
			this->weights.push_back ( Gw_3_B );
			break;                                                                             }
		case tag::tri_4 :  // Gauss quadrature with four points on a triangle
		{ FiniteElement::WithMaster::Triangle * tri_fe =
				dynamic_cast < FiniteElement::WithMaster::Triangle * > ( fe.core );
			assert ( tri_fe );
			this->points.reserve(4);  this->weights.reserve(4);
			// points
			Cell Gauss_4_c ( tag::vertex );  xi(Gauss_4_c) = 1./3.;  eta(Gauss_4_c) = 1./3.;
			Cell Gauss_4_O ( tag::vertex );  xi(Gauss_4_O) = 0.2;    eta(Gauss_4_O) = 0.2;
			Cell Gauss_4_A ( tag::vertex );  xi(Gauss_4_A) = 0.6;    eta(Gauss_4_A) = 0.2;
			Cell Gauss_4_B ( tag::vertex );  xi(Gauss_4_B) = 0.2;    eta(Gauss_4_B) = 0.6;
			// weights
			double Gw_4_c = -27./96.;  // yes, negative weight
			double Gw_4_O = 25./96., Gw_4_A = 25./96., Gw_4_B = 25./96.;
			// in the book of Oden we have the same weights but a different distribution of points
			// instead of  0.2 0.2 0.6  we have  2/15 2/15 11/15
			// maybe the points are not unique ?
			this->points.push_back ( Gauss_4_c );
			this->weights.push_back ( Gw_4_c );
			this->points.push_back ( Gauss_4_O );
			this->weights.push_back ( Gw_4_O );
			this->points.push_back ( Gauss_4_A );
			this->weights.push_back ( Gw_4_A );
			this->points.push_back ( Gauss_4_B );
			this->weights.push_back ( Gw_4_B );
			break;                                                                                }
		case tag::tri_4_Oden :
		// Gauss quadrature with four points on a triangle, points presented in the book
		// E.B. Becker, G.F. Carey, J.T. Oden, Finite Elements, an introduction, vol 1
		{ FiniteElement::WithMaster::Triangle * tri_fe =
				dynamic_cast < FiniteElement::WithMaster::Triangle * > ( fe.core );
			assert ( tri_fe );
			this->points.reserve(4);  this->weights.reserve(4);
			// points
			Cell Gauss_4_c ( tag::vertex );  xi(Gauss_4_c) =  1./ 3.;  eta(Gauss_4_c) =  1./ 3.;
			Cell Gauss_4_O ( tag::vertex );  xi(Gauss_4_O) =  2./15.;  eta(Gauss_4_O) =  2./15.;
			Cell Gauss_4_A ( tag::vertex );  xi(Gauss_4_A) = 11./15.;  eta(Gauss_4_A) =  2./15.;
			Cell Gauss_4_B ( tag::vertex );  xi(Gauss_4_B) =  2./15.;  eta(Gauss_4_B) = 11./15.;
			// weights
			double Gw_4_c = -27./96.;  // yes, negative weight
			double Gw_4_O = 25./96., Gw_4_A = 25./96., Gw_4_B = 25./96.;
			this->points.push_back ( Gauss_4_c );
			this->weights.push_back ( Gw_4_c );
			this->points.push_back ( Gauss_4_O );
			this->weights.push_back ( Gw_4_O );
			this->points.push_back ( Gauss_4_A );
			this->weights.push_back ( Gw_4_A );
			this->points.push_back ( Gauss_4_B );
			this->weights.push_back ( Gw_4_B );
			break;                                                                                }
		case tag::tri_6 :  // Gauss quadrature with six points on a triangle
		{	FiniteElement::WithMaster::Triangle * tri_fe =
				dynamic_cast < FiniteElement::WithMaster::Triangle * > ( fe.core );
			assert ( tri_fe );
			this->points.reserve(6);  this->weights.reserve(6);
			// points
			Cell Gauss_6_O ( tag::vertex );   xi (Gauss_6_O)  = 0.091576213509771;
			                                  eta(Gauss_6_O)  = 0.091576213509771;
			Cell Gauss_6_A ( tag::vertex );   xi (Gauss_6_A)  = 0.816847572980459;
			                                  eta(Gauss_6_A)  = 0.091576213509771;
			Cell Gauss_6_B ( tag::vertex );   xi (Gauss_6_B)  = 0.091576213509771;
			                                  eta(Gauss_6_B)  = 0.816847572980459;
			Cell Gauss_6_OA ( tag::vertex );  xi (Gauss_6_OA) = 0.445948490915965;
			                                  eta(Gauss_6_OA) = 0.108103018168070;
			Cell Gauss_6_OB ( tag::vertex );  xi (Gauss_6_OB) = 0.108103018168070;
			                                  eta(Gauss_6_OB) = 0.445948490915965;
			Cell Gauss_6_AB ( tag::vertex );  xi (Gauss_6_AB) = 0.445948490915965;
			                                  eta(Gauss_6_AB) = 0.445948490915965;
			// weights
			double Gw_6_O = 0.054975871827661, Gw_6_A = 0.054975871827661, Gw_6_B = 0.054975871827661;
			double Gw_6_OA = 0.1116907948390055, Gw_6_OB = 0.1116907948390055, Gw_6_AB = 0.1116907948390055;
			this->points.push_back ( Gauss_6_O );
			this->weights.push_back ( Gw_6_O );
			this->points.push_back ( Gauss_6_A );
			this->weights.push_back ( Gw_6_A );
			this->points.push_back ( Gauss_6_B );
			this->weights.push_back ( Gw_6_B );
			this->points.push_back ( Gauss_6_OA );
			this->weights.push_back ( Gw_6_OA );
			this->points.push_back ( Gauss_6_OB );
			this->weights.push_back ( Gw_6_OB );
			this->points.push_back ( Gauss_6_AB );
			this->weights.push_back ( Gw_6_AB );
			break;                                                                                           }
		case tag::quad_4 :  // Gauss quadrature with four points on a quadrangle
		{ FiniteElement::WithMaster::Quadrangle * quadr_fe =
				dynamic_cast < FiniteElement::WithMaster::Quadrangle * > ( fe.core );
			assert ( quadr_fe );
			this->points.reserve(4);  this->weights.reserve(4);
			const double sqrt_of_one_third = sqrt(1./3.);
			Cell Gauss_4_a ( tag::vertex );   xi (Gauss_4_a)  = -sqrt_of_one_third;
			                                  eta(Gauss_4_a)  = -sqrt_of_one_third;
			Cell Gauss_4_b ( tag::vertex );   xi (Gauss_4_b)  = -sqrt_of_one_third;
			                                  eta(Gauss_4_b)  =  sqrt_of_one_third;
			Cell Gauss_4_c ( tag::vertex );   xi (Gauss_4_c)  =  sqrt_of_one_third;
			                                  eta(Gauss_4_c)  = -sqrt_of_one_third;
			Cell Gauss_4_d ( tag::vertex );   xi (Gauss_4_d)  =  sqrt_of_one_third;
			                                  eta(Gauss_4_d)  =  sqrt_of_one_third;
			this->points.push_back ( Gauss_4_a );
			this->weights.push_back ( 1. );
			this->points.push_back ( Gauss_4_b );
			this->weights.push_back ( 1. );
			this->points.push_back ( Gauss_4_c );
			this->weights.push_back ( 1. );
			this->points.push_back ( Gauss_4_d );
			this->weights.push_back ( 1. );
			break;                                                                    }
		case tag::quad_9 :  // Gauss quadrature with nine points on a quadrangle
		{	FiniteElement::WithMaster::Quadrangle * quadr_fe =
				dynamic_cast < FiniteElement::WithMaster::Quadrangle * > ( fe.core );
			assert ( quadr_fe );
			this->points.reserve(9);  this->weights.reserve(9);
			// points :
			const double sqrt3over5 = sqrt(0.6);
			Cell Gauss_9_a ( tag::vertex );   xi (Gauss_9_a)  = -sqrt3over5;
			                                  eta(Gauss_9_a)  = -sqrt3over5;
			Cell Gauss_9_b ( tag::vertex );   xi (Gauss_9_b)  =  0.;
			                                  eta(Gauss_9_b)  = -sqrt3over5;
			Cell Gauss_9_c ( tag::vertex );   xi (Gauss_9_c)  =  sqrt3over5;
			                                  eta(Gauss_9_c)  = -sqrt3over5;
			Cell Gauss_9_d ( tag::vertex );   xi (Gauss_9_d)  = -sqrt3over5;
			                                  eta(Gauss_9_d)  =  0.;
			Cell Gauss_9_e ( tag::vertex );   xi (Gauss_9_e)  =  0.;
			                                  eta(Gauss_9_e)  =  0.;
			Cell Gauss_9_f ( tag::vertex );   xi (Gauss_9_f)  =  sqrt3over5;
			                                  eta(Gauss_9_f)  =  0.;
			Cell Gauss_9_g ( tag::vertex );   xi (Gauss_9_g)  = -sqrt3over5;
			                                  eta(Gauss_9_g)  =  sqrt3over5;
			Cell Gauss_9_h ( tag::vertex );   xi (Gauss_9_h)  =  0.;
			                                  eta(Gauss_9_h)  =  sqrt3over5;
			Cell Gauss_9_i ( tag::vertex );   xi (Gauss_9_i)  =  sqrt3over5;
			                                  eta(Gauss_9_i)  =  sqrt3over5;
			//  a{-sqrt3over5,-sqrt3over5}, b{0, -sqrt3over5}, c{sqrt3over5, -sqrt3over5},
			//  d{-sqrt3over5,0}, e{0,0}, f{sqrt3over5, 0}, g{-sqrt3over5,sqrt3over5},
			//  h{0,sqrt3over5}, i{sqrt3over5, sqrt3over5}
			// weights :
			//  25./81., 40./81., 25./81., 40./81, 64./81., 40./81, 25./81., 40./81, 25./81.
			this->points.push_back ( Gauss_9_a );
			this->weights.push_back ( 25./81. );
			this->points.push_back ( Gauss_9_b );
			this->weights.push_back ( 40./81. );
			this->points.push_back ( Gauss_9_c );
			this->weights.push_back ( 25./81. );
			this->points.push_back ( Gauss_9_d );
			this->weights.push_back ( 40./81. );
			this->points.push_back ( Gauss_9_e );
			this->weights.push_back ( 64./81. );
			this->points.push_back ( Gauss_9_f );
			this->weights.push_back ( 40./81. );
			this->points.push_back ( Gauss_9_g );
			this->weights.push_back ( 25./81. );
			this->points.push_back ( Gauss_9_h );
			this->weights.push_back ( 40./81. );
			this->points.push_back ( Gauss_9_i );
			this->weights.push_back ( 25./81. );
			break;                                                                   }
		default :	assert ( false );
	}  // end of  switch ( q )
}
	
//-----------------------------------------------------------------------------------------//

double Integrator::Gauss::action ( const Function & f, const FiniteElement & fe )
// virtual from Integrator::Core

// assumes the finite element is already docked on a cell
// thus, fe_core->transf is well defined

{	FiniteElement::WithMaster * fe_core = dynamic_cast < FiniteElement::WithMaster * > ( fe.core );
	assert ( fe_core );
	Function::Diffeomorphism * tran = Function::core_to_diffeom ( fe_core->transf.core );
	
	// Function ff = f.replace ( tran->geom_coords, tran->back_geom_coords );
	// nope. we assume this has already been done

	Function ff = f;

	double res = 0.;
	std::vector<double>::iterator it_weight = this->weights.begin();
	for ( std::vector<Cell>::iterator it_point = this->points.begin();
	      it_point != this->points.end(); it_point++                   )
	{	assert ( it_weight != this->weights.end() );
		Cell Gauss_point = *it_point;
		double w = *it_weight;
		res += w * ff(Gauss_point) * tran->det(Gauss_point);
		it_weight++;                                   }
	assert ( it_weight == this->weights.end() );
	return res;                                                            }

//-----------------------------------------------------------------------------------------//

void FiniteElement::WithMaster::Triangle::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core, through FiniteElement::WithMaster

{	assert ( cll.dim() == 2 );
	this->docked_on = cll;
	CellIterator it = cll.boundary().iter_over ( tag::vertices );
	it.reset();  assert ( it.in_range() );  Cell P = *it;
	it++;  assert ( it.in_range() );  Cell Q = *it;
	it++;  assert ( it.in_range() );  Cell R = *it;
	it++;  assert ( not it.in_range() );
	
	Function xy = Manifold::working.coordinates();
	assert ( xy.nb_of_components() == 2 );
	Function x = xy[0],  y = xy[1];
	double xP = x(P), xQ = x(Q), xR = x(R);
	double yP = y(P), yQ = y(Q), yR = y(R);

	Function xi_eta = this->master_manif.coordinates();
	assert ( xi_eta.nb_of_components() == 2 );
	Function xi = xi_eta[0], eta = xi_eta[1];

	Function one_m_xi_m_eta = 1.-xi-eta;
	Function x_c = xP * one_m_xi_m_eta + xQ * xi + xR * eta;
	Function y_c = yP * one_m_xi_m_eta + yQ * xi + yR * eta;
	
	
	// xi, eta and 1-xi-eta are a base of functions defined on the master element
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix

	this->transf = Function ( tag::diffeomorphism, xy, xi_eta, x_c && y_c );

	this->base_fun_1.clear();
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( one_m_xi_m_eta, tag::composed_with, this->transf ) ) );
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( xi, tag::composed_with, this->transf ) ) );
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( eta, tag::composed_with, this->transf ) ) );              }
	
//-----------------------------------------------------------------------------------------//

void FiniteElement::WithMaster::Quadrangle::dock_on ( const Cell & cll )
// virtual from FiniteElement::Core, through FiniteElement::WithMaster

{	assert ( cll.dim() == 2 );
	this->docked_on = cll;
	CellIterator it = cll.boundary().iter_over ( tag::vertices );
	it.reset();  assert ( it.in_range() );  Cell P = *it;
	it++;  assert ( it.in_range() );  Cell Q = *it;
	it++;  assert ( it.in_range() );  Cell R = *it;
	it++;  assert ( it.in_range() );  Cell S = *it;
	it++;  assert ( not it.in_range() );
	
	Function xy = Manifold::working.coordinates();
	assert ( xy.nb_of_components() == 2 );
	Function x = xy[0],  y = xy[1];
	double xP = x(P), xQ = x(Q), xR = x(R), xS = x(S);
	double yP = y(P), yQ = y(Q), yR = y(R), yS = y(S);

	Function xi_eta = this->master_manif.coordinates();
	assert ( xi_eta.nb_of_components() == 2 );
	Function xi = xi_eta[0], eta = xi_eta[1];
	
	Function psiP = (1.-xi)*(1.-eta), psiQ = (1.+xi)*(1.-eta),
		       psiR = (1.+xi)*(1.+eta), psiS = (1.-xi)*(1.+eta);
	// psi* are a base of functions defined on the master element
	// each takes value 4 in the associated vertex, zero on the others
	// we may need to differentiate them with respect to x and y (physical coordinates)
	// which is equivalent to differentiating x_c and y_c with respect to xi and eta
	// and then taking the inverse matrix
	
	Function x_c = ( xP * psiP + xQ * psiQ + xR * psiR + xS * psiS ) / 4.;
	Function y_c = ( yP * psiP + yQ * psiQ + yR * psiR + yS * psiS ) / 4.;

	this->transf = Function ( tag::diffeomorphism, xy, xi_eta, x_c && y_c );

	this->base_fun_1.clear();
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( P.core, Function ( psiP/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( Q.core, Function ( psiQ/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( R.core, Function ( psiR/4., tag::composed_with, this->transf ) ) );
	this->base_fun_1.insert ( std::pair < Cell::Core*, Function >
	       ( S.core, Function ( psiS/4., tag::composed_with, this->transf ) ) );     }

//-----------------------------------------------------------------------------------------//

