
// maniFEM manifold.cpp 2019.11.05

#include "manifold.h"

using namespace maniFEM;

	
Manifold Manifold::current ( tag::ghost );
// anything would do, the user must set this variable before anything else
// by simply declaring any Manifold object (constructor assigns to Manifold::current)

	
Function Manifold::Euclid::build_coord_func  // maybe change to not virtual ?
( const tag::lagrange &, const tag::OfDegree &, size_t deg )
// virtual from Function::Core
	
{	assert ( deg == 1 );
	Field::Core * field;
	Function::Core * func;
	if ( this->dim == 1 )
	{	Field::Scalar * field_scalar = new Field::Scalar
			( tag::lives_on_positive_cells, tag::of_dim, 0 );
		field = field_scalar;
		func = new Function::CoupledWithField::Scalar ( field_scalar );  }
	else
	{	assert ( this->dim > 1 );
		Field::Block * field_block = new Field::Block
			( tag::lives_on_positive_cells, tag::of_dim, 0, tag::has_size, this->dim );
		field = field_block;
	  func = new Function::CoupledWithField::Vector ( field_block );                }
	this->coord_field = field;
	this->coord_func = Function ( tag::whose_core_is,
                                std::shared_ptr < Function::Core > ( func ) );
	return this->coord_func;                                                            }


Function Manifold::Implicit::build_coord_func
( const tag::lagrange &, const tag::OfDegree &, size_t deg )
// virtual from Function::Core
{	assert ( false );  }

	
Function Manifold::Euclid::get_coord_func ( ) const  // virtual from Function::Core
{	return this->coord_func;  }
	
Function Manifold::Implicit::get_coord_func ( ) const  // virtual from Function::Core
{	return this->surrounding_space.coordinates();  }
	


// P = sA + sB,  s+t == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B);  }	

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster

void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const
//  virtual from Manifold::Core

{	Function coord = this->get_coord_func();
	std::shared_ptr < Function::Scalar > coord_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B) );
		return;                                            } 
	std::shared_ptr < Function::Vector > coord_vector =
		std::dynamic_pointer_cast < Function::Vector > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	std::shared_ptr < Function::Scalar > coord_i = coord_vector->component(i);
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) );                                   }  }


// P = sA + sB + uC + vD,  s+t+u+v == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B) +
		              u * coord[i](C) + v * coord[i](D);   }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster

void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const
//  virtual from Manifold::Core

{	Function coord = this->get_coord_func();
	std::shared_ptr < Function::Scalar > coord_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B) +
			     u * coord_scalar->get_value_on_cell(C) +
			     v * coord_scalar->get_value_on_cell(D)   );
		return;                                                        } 
	std::shared_ptr < Function::Vector > coord_vector =
		std::dynamic_pointer_cast < Function::Vector > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	std::shared_ptr < Function::Scalar > coord_i = coord_vector->component(i);
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) +
			     u * coord_i->get_value_on_cell(C) +
			     v * coord_i->get_value_on_cell(D)   );                                 }   }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1

void Manifold::Euclid::pretty_interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D,
                  double w, const Cell & E, double z, const Cell & F ) const

{	Function coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
		coord[i](P) = s * coord[i](A) + t * coord[i](B) +
		              u * coord[i](C) + v * coord[i](D) +
		              w * coord[i](E) + z * coord[i](F);   }

// code below does the same as code above
// above is pretty, slightly slower - below is ugly, slightly faster

void Manifold::Euclid::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const
//  virtual from Manifold::Core

{	Function coord = this->get_coord_func();
	std::shared_ptr < Function::Scalar > coord_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( coord.core );
	if ( coord_scalar )
	{	assert ( coord.nb_of_components() == 1 );
	  coord_scalar->set_value_on_cell
			( P, s * coord_scalar->get_value_on_cell(A) +
			     t * coord_scalar->get_value_on_cell(B) +
			     u * coord_scalar->get_value_on_cell(C) +
			     v * coord_scalar->get_value_on_cell(D) +
			     w * coord_scalar->get_value_on_cell(E) +
			     z * coord_scalar->get_value_on_cell(F)   );
		return;                                             } 
	std::shared_ptr < Function::Vector > coord_vector =
		std::dynamic_pointer_cast < Function::Vector > ( coord.core );
	assert ( coord_vector );
	size_t n = coord.nb_of_components();
	for ( size_t i = 0; i < n; i++ )
	{	std::shared_ptr < Function::Scalar > coord_i = coord_vector->component(i);
		coord_i->set_value_on_cell
			( P, s * coord_i->get_value_on_cell(A) +
			     t * coord_i->get_value_on_cell(B) +
				   u * coord_i->get_value_on_cell(C) +
			     v * coord_i->get_value_on_cell(D) +
			     w * coord_i->get_value_on_cell(E) +
			     z * coord_i->get_value_on_cell(F)   );                                  }   }


// P = sA + sB,  s+t == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B );
	this->project ( P );                                          }


// P = sA + sB + uC + vD,  s+t == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D );
	this->project ( P );                                                      }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
void Manifold::Implicit::interpolate ( Cell::Positive::Vertex * P,
  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const
//  virtual from Manifold::Core

{	this->surrounding_space.core->interpolate ( P, s, A, t, B, u, C, v, D, w, E, z, F );
	this->project ( P );                                                                  }


void Manifold::Euclid::project ( Cell::Positive::Vertex * P_c ) const
{	}
	
void Manifold::Implicit::project ( Cell::Positive::Vertex * P_c ) const

// just a few steps of Newton's method for under-determined (systems of) equations

{	const Function & coord = this->get_coord_func();
	size_t n = coord.nb_of_components();
	const Function & lev_func = this->level_function;
	assert ( lev_func.nb_of_components() == 1 );  // for now, just one equation
	const Function & grad_lev = this->grad_lev_func;
	assert ( grad_lev.nb_of_components() == n );
	const Cell P ( tag::whose_core_is, P_c );  // temporary wrapper for P_c
	for ( short int k = 0; k < Manifold::Implicit::steps_for_Newton; k++ )
	// we move doubles around a lot
	// how to do it faster ?
	// somehow bind references to coord_at_P to 'coord'
	{	std::vector<double> coord_at_P = coord ( P );
		double lev_at_P = lev_func ( P );
		std::vector<double> grad_lev_at_P = grad_lev ( P );
		double norm2 = 0.;
		for ( size_t i = 0; i < n; i++ )
			norm2 += grad_lev_at_P[i]*grad_lev_at_P[i];
		double coef = lev_at_P / norm2;
		for ( size_t i = 0; i < n; i++ )
			coord_at_P[i] -= coef * grad_lev_at_P[i];
		coord ( P ) = coord_at_P;                               }                    }
