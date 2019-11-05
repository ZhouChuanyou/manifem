
// maniFEM manifold.h 2019.11.05

#ifndef MANIFEM_MANIFOLD_H
#define MANIFEM_MANIFOLD_H

#include <iostream>
#include <vector>
#include <forward_list>
#include "assert.h"

#include "mesh.h"
#include "field.h"
#include "function.h"


namespace maniFEM {

namespace tag {

	struct euclid { };  static const euclid Euclid;
	struct Implicit { };  static const Implicit implicit;
}


class Manifold

{	public :

	class Core;
	
	Manifold::Core * core;

	inline Manifold ( const tag::WhoseCoreIs &, Manifold::Core * c )
	:	core ( c )
	{	assert ( c );
		Manifold::current = *this;  }

	inline Manifold ( const Manifold & m ) : core ( m.core ) { }

	inline Manifold ( const Manifold && m ) : core ( m.core ) { }
	// is this the right way to do it ?

	inline Manifold ( const tag::Ghost & ) : core ( nullptr ) { }

	inline Manifold ( const tag::euclid &, const tag::OfDimension &, size_t dim );

	inline Manifold ( const tag::Implicit &, const Manifold &, const Function & );
	
	inline Manifold & operator= ( const Manifold & m )
	{	core = m.core;
		return *this;   }
	
	inline Manifold & operator= ( const Manifold && m )
	{	core = m.core;  // is this the right way to do it ?
		return *this;   }
	
	inline Function coordinates ( ) const;

	inline Function build_coordinate_system
	( const tag::lagrange &, const tag::OfDegree &, size_t d );

	// only a ghost manifold has null core
	inline bool exists ( ) const { return core != nullptr; } 
	
	// P = sA + sB,  s+t == 1
	inline void interpolate
	( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	inline void interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
	  double w, const Cell & E, double z, const Cell & F ) const;

	inline void project ( const Cell & ) const;
	
	inline Manifold implicit ( const Function::Equation eq ) const;

	class Euclid;  class Implicit;

	static Manifold current;

};

class Manifold::Core

{	public :

	size_t dim;

	// P = sA + sB,  s+t == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const = 0;

	// P = sA + sB + uC + vD,  s+t+u+v == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D ) const = 0;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
	virtual void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const = 0;

	virtual Function build_coord_func
		( const tag::lagrange &, const tag::OfDegree &, size_t d ) = 0;
	
	virtual Function get_coord_func ( ) const = 0;
	
	virtual void project ( Cell::Positive::Vertex * ) const = 0;
};


inline Function Manifold::coordinates ( ) const
{	assert ( this->core );
	return this->core->get_coord_func();  }

inline Function Manifold::build_coordinate_system
( const tag::lagrange &, const tag::OfDegree &, size_t d )
{	return core->build_coord_func ( tag::Lagrange, tag::of_degree, d );  }


// P = sA + sB,  s+t == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const
{	assert ( P.dim() == 0 );  assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( P.is_positive() );  assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
	  s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core );  }


// P = sA + sB + uC + vD,  s+t+u+v == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D ) const
{	assert ( P.dim() == 0 );
	assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( C.dim() == 0 );  assert ( D.dim() == 0 );
	assert ( P.is_positive() );
	assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( C.is_positive() );  assert ( D.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
	  s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core,
	  u, ( Cell::Positive::Vertex * ) C.core, v, ( Cell::Positive::Vertex * ) D.core );  }


// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1
inline void Manifold::interpolate
( const Cell & P, double s, const Cell & A, double t, const Cell & B,
                  double u, const Cell & C, double v, const Cell & D,
                  double w, const Cell & E, double z, const Cell & F ) const
{	assert ( P.dim() == 0 );
	assert ( A.dim() == 0 );  assert ( B.dim() == 0 );
	assert ( C.dim() == 0 );  assert ( D.dim() == 0 );
	assert ( E.dim() == 0 );  assert ( F.dim() == 0 );
	assert ( P.is_positive() );
	assert ( A.is_positive() );  assert ( B.is_positive() );
	assert ( C.is_positive() );  assert ( D.is_positive() );
	assert ( E.is_positive() );  assert ( F.is_positive() );
	assert ( s >= 0. );  assert ( s <= 1. );
	assert ( t >= 0. );  assert ( t <= 1. );
	assert ( u >= 0. );  assert ( u <= 1. );
	assert ( v >= 0. );  assert ( v <= 1. );
	assert ( w >= 0. );  assert ( w <= 1. );
	assert ( z >= 0. );  assert ( z <= 1. );
	// cannot assert the sum is 1. due to round-off errors
	this->core->interpolate ( ( Cell::Positive::Vertex * ) P.core,
	  s, ( Cell::Positive::Vertex * ) A.core, t, ( Cell::Positive::Vertex * ) B.core,
	  u, ( Cell::Positive::Vertex * ) C.core, v, ( Cell::Positive::Vertex * ) D.core,
	  w, ( Cell::Positive::Vertex * ) E.core, z, ( Cell::Positive::Vertex * ) F.core );  }


inline void Manifold::project ( const Cell & cll ) const
{	assert ( cll.is_positive () );
	assert ( cll.dim() == 0 );
	this->core->project ( ( Cell::Positive::Vertex * ) cll.core );  }


class Manifold::Euclid : public Manifold::Core

{	public :

	size_t dim;
	
	Field::Core * coord_field { nullptr };
	Function coord_func { tag::non_existent };

	inline Euclid ( size_t d );

	// P = sA + sB,  s+t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void pretty_interpolate
		( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;
	void pretty_interpolate
	(	const Cell & P, double s, const Cell & A, double t, const Cell & B,
		                double u, const Cell & C, double v, const Cell & D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;
	void pretty_interpolate ( const Cell & P, double s, const Cell & A,
	  double t, const Cell & B, double u, const Cell & C, double v, const Cell & D,
	  double w, const Cell & E, double z, const Cell & F ) const;

  Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d ) ;
  // virtual from Function::Core
	
  Function get_coord_func ( ) const;  // virtual from Function::Core
	
	void project ( Cell::Positive::Vertex * ) const;  // virtual from Manifold::Core

};  // end of class Manifold::Euclid


// when we first declare a Manifold::Euclid, coord_field and coord_func will be nullptr
// immediately after declaring the manifold (in main), we must declare the coordinates
// to be of type Lagrange degree 1 for instance
//    	Manifold RR2 ( tag::Euclid, 2);
//      Function xy = RR2.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );

inline Manifold::Manifold ( const tag::euclid &, const tag::OfDimension &, size_t d )
:	Manifold ( tag::whose_core_is, new Manifold::Euclid ( d ) )
{	assert ( d > 0 );  }


inline Manifold::Euclid::Euclid ( size_t d )
:	Manifold::Core(), dim { d }
{	assert ( d > 0 );  }


class Manifold::Implicit : public Manifold::Core

// a submanifold of a Manifold::Euclid defined by one or more equations

{	public :

	Manifold surrounding_space;
	Function level_function, grad_lev_func;

	// the projection will be done by means of the Newton method
	static const short int steps_for_Newton = 10;
	
	inline Implicit ( const Manifold & s, const Function & f );
	
	// P = sA + sB,  s+t == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B ) const;
	void pretty_interpolate
		( const Cell & P, double s, const Cell & A, double t, const Cell & B ) const;

	// P = sA + sB + uC + vD,  s+t+u+v == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D  ) const;

	// P = sA + sB + uC + vD + wE + zF,  s+t+u+v+w+z == 1     virtual from Manifold::Core
	void interpolate ( Cell::Positive::Vertex * P,
	  double s, Cell::Positive::Vertex * A, double t, Cell::Positive::Vertex * B,
	  double u, Cell::Positive::Vertex * C, double v, Cell::Positive::Vertex * D,
	  double w, Cell::Positive::Vertex * E, double z, Cell::Positive::Vertex * F ) const;

	Function build_coord_func ( const tag::lagrange &, const tag::OfDegree &, size_t d );
	//   virtual from Function::Core
	
	Function get_coord_func ( ) const;  // virtual from Function::Core

	void project ( Cell::Positive::Vertex * ) const;  // virtual from Manifold::Core
	
};  // end of class Manifold::Implicit


inline Manifold::Manifold ( const tag::Implicit &, const Manifold & m, const Function & f )
:	Manifold ( tag::whose_core_is, new Manifold::Implicit ( m, f ) )
{	}


Manifold::Implicit::Implicit ( const Manifold & m, const Function & f )
:	surrounding_space ( m ), level_function ( f ),
	grad_lev_func ( tag::non_existent )  // temporarily empty gradient

{	Function coord = m.coordinates();
	size_t n = coord.nb_of_components();
	Function::Aggregate * grad = new Function::Aggregate ( tag::has_size, n );
	for ( size_t i = 0; i < n; i++ )
	{	Function partial = f.deriv(coord[i]);
		std::shared_ptr<Function::Scalar> partial_c =
			std::dynamic_pointer_cast<Function::Scalar> ( partial.core );
		assert ( partial_c );
		grad->components[i] = partial_c;                                }
	this->grad_lev_func.core = std::shared_ptr<Function::Core> ( grad );        }


inline Manifold Manifold::implicit ( const Function::Equation eq ) const
{	return Manifold ( tag::implicit, *this, eq.lhs );  }


inline void Cell::project () const
{	this->project ( tag::onto, Manifold::current );  }
	
inline void Cell::project ( const tag::Onto &, const Manifold m ) const
{	Manifold::Implicit * impl = dynamic_cast<Manifold::Implicit*> ( m.core );
	assert ( impl );
	Cell::Positive::Vertex * cll = dynamic_cast<Cell::Positive::Vertex*> ( this->core );
	assert ( cll );
	impl->project ( cll );                                                                }
		

}  // namespace maniFEM

#endif
// ifndef MANIFEM_MANIFOLD_H




