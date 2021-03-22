
// function.h 2021.03.21

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

#ifndef MANIFEM_FUNCTION_H
#define MANIFEM_FUNCTION_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <memory>
#include "math.h"
#include "assert.h"

#include "mesh.h"
#include "field.h"


namespace maniFEM {


namespace tag
{	struct Diffeomorphism { };  static const Diffeomorphism diffeomorphism;
	struct ComposedWith { };  static const ComposedWith composed_with;
	struct Iff { };  static const Iff iff;
	struct LessThan { };  static const LessThan less_than;
	struct IfLessThan { };  static const IfLessThan if_less_than;
	struct Otherwise { };  static const Otherwise otherwise;                  }
	

class Manifold;
	

class Function

// just a a thin wrapper around Function::Core
// I guess this is what they call "delegation" ...
// useful for arithmetic operators and for vector functions

{	public :

	class Core;

	#ifndef NDEBUG
	static std::map < const Function::Core*, std::string > name;
	#endif

	std::shared_ptr < Function::Core > core;

	inline Function ( const tag::NonExistent & ) : core { nullptr } { };

	inline Function ( const tag::WhoseCoreIs &, Function::Core * c );

	inline Function ( const Function & f )
	:	Function ( tag::whose_core_is, f.core )  { }

	inline Function ( const tag::HasSize &, size_t s );

	inline Function ( double c );

	inline Function ( const tag::Diffeomorphism &, const Function & geom_coords,
                    const Function & master_coords, const Function & geom_back_coords );

	inline Function ( const Function & expr,
                    const tag::ComposedWith &, const Function & diffeo );
	// an 'expr'ession involving master coordinates ( e.g.  1. - xi - eta )
	// composed with a 'diffeo'morphism sending it in the physical space

	inline ~Function();

	inline Function & operator= ( const Function & m );
	
	inline size_t nb_of_components ( ) const;

	inline Function operator[] ( size_t ) const;

	class TakenOnCell;
	inline Function::TakenOnCell operator() ( const Cell & cll ) const;

	static Function Zero, One, MinusOne;

	inline Function deriv ( const Function & x ) const;  // derivative with respect to x

	inline Function replace ( const Function & x, const Function & y ) const;
	// in an expression, replace x by y

	struct Catalogue
	{	// std::shared_ptr < Function::Core > self;
		Function * self;
		std::map < Function::Core *, std::shared_ptr < Function::Core > > deriv;  };
	
	#ifndef NDEBUG
	inline std::string repr ( ) const;
	enum From { from_void, from_sum, from_product, from_power, from_function };
	#endif

	class Scalar;  class ArithmeticExpression;  class Constant;
	class Vector;  class Aggregate;  class CoupledWithField;
	class Sum;  class Product;  class Power;  class Sin;  class Cos;  class Step;
	class Diffeomorphism;  class Composition;
	class Equality;

	static inline Function::Scalar * core_to_scalar ( Function::Core * f );
	static inline Function::Vector * core_to_vector ( Function::Core * f );
	static inline Function::Constant * core_to_constant ( Function::Core * f );
	static inline Function::Sum * core_to_sum ( Function::Core * f );
	static inline Function::Product * core_to_product ( Function::Core * f );
	static inline Function::Diffeomorphism * core_to_diffeom ( Function::Core * f );

};  // end of  class Function


//-----------------------------------------------------------------------------------------//

class Function::Core

// just a base for classes like Constant, Sum, Product and many others
// we use dynamic_cast

{	public :

	// Manifold * manifold;

	size_t number_of_wrappers { 0 };
	
	// inline Core ( ) : number_of_wrappers { 0 } { };

	virtual ~Core ( ) { };

	inline void release ( );

	virtual size_t nb_of_components ( ) const = 0;

	virtual Function component ( size_t i ) = 0;

	virtual Function deriv ( Function ) const = 0;

	virtual Function replace ( const Function & x, const Function & y ) = 0;
	// in an expression, replace x by y

	#ifndef NDEBUG	
	virtual std::string repr ( const Function::From & from = Function::from_void ) const = 0;
	#endif
};

//-----------------------------------------------------------------------------------------//

inline void Function::Core::release ( )

{	assert ( this->number_of_wrappers > 0 );
	// std::cout << "Function destructor " << this->number_of_wrappers << " "
	//					<< this << " " << Function::name[this] << std::endl;
	this->number_of_wrappers --;
	if ( this->number_of_wrappers == 0 ) delete this;  }


inline Function::Function ( const tag::WhoseCoreIs &, Function::Core * c )
:	core ( c )
{	assert ( c );  c->number_of_wrappers ++;  }
	// std::cout << "Function constructor " << c->number_of_wrappers << " " << c << std::endl;  }

inline Function::~Function ()
{	if ( this->core ) this->core->release();  }

inline Function & Function::operator= ( const Function & m )
{	if ( this->core ) this->core->release();
	this->core = m.core;
	this->core->number_of_wrappers ++;
	return *this;                             }
	
inline Function Function::replace ( const Function & x, const Function & y ) const
{	return this->core->replace ( x, y );  }

//-----------------------------------------------------------------------------------------//

class Function::Scalar : public Function::Core
	
{	public :

	size_t nb_of_components ( ) const;  // virtual from Function::Core, here returns 1

	Function component ( size_t i ); // virtual from Function::Core

	virtual double set_value_on_cell ( Cell::Core *, const double & ) = 0;
	// assign a numeric value to the function on the cell and return that value
	
	virtual double get_value_on_cell ( Cell::Core * ) const = 0;

	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//    stay pure virtual from Function::Core
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core
	#endif
};


//-----------------------------------------------------------------------------------------//

class Function::ArithmeticExpression : public Function::Scalar
	
// base for classes like Constant, Sum, Product
// any function for which 'set_value_on_cell' does not make sense

{	public :

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	// here execution forbidden
	
	// double get_value_on_cell ( Cell::Core * )  stays pure virtual from Function::Scalar

	// Function deriv ( Function )
	// Function replace ( const Function & x, const Function & y )
	//   stay pure virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
	//   stays pure virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Constant : public Function::ArithmeticExpression
	
{	public :

	double val;

	inline Constant ( double c )
	:	val { c } { }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function::Function ( double c )
:	Function ( tag::whose_core_is, new Function::Constant ( c ) )
{	}

//-----------------------------------------------------------------------------------------//

class Function::Sum : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < Function > terms;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Product : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < Function > factors;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Power : public Function::ArithmeticExpression
	
{	public :

	Function base;
	double exponent;

	inline Power ( Function b, double e )
	:	base { b }, exponent { e }
	{	assert ( dynamic_cast < Function::Scalar * > ( b.core ) );  }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Sin : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer

	inline Sin ( const Function & b ) : base { b }  { }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function sin ( const Function & f )
{	return Function ( tag::whose_core_is,
		new Function::Sin ( f ) );  }

//-----------------------------------------------------------------------------------------//

class Function::Cos : public Function::ArithmeticExpression
	
{	public :

	Function base;  // we use wrapper as a pointer

	inline Cos ( const Function & b ) : base { b }  { }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function cos ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Cos ( f ) );  }

//-----------------------------------------------------------------------------------------//

class Function::Step : public Function::ArithmeticExpression

// piecewise constant functions 

{	public :

	Function arg;  // we use wrappers as pointers
	// 'arg' is used for testing inequalities

	std::vector < Function > values;
	std::vector < double > cuts;  //  values.size() == cuts.size() + 1
	// 'cuts' contains values ordered increasingly
	// for instance, suppose values=={x,1/x,3} and cuts=={-1,1}
	// this means this(x)=x for x<-1, this(x)=1/x for -1<x<1, this(x)=3 for x>1

	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c,
	              const Function & v2, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2 }, cuts { c }
	{	}
	
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
	              const Function & v3, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3 }, cuts { c1, c2 }
	{	assert ( c1 < c2 );  }
	
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
                const Function & v3, const tag::IfLessThan, double c3,
	              const Function & v4, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3, v4 }, cuts { c1, c2, c3 }
	{	assert ( c1 < c2 );  assert ( c2 < c3 );  }
	
	inline Step ( const Function & v1, const tag::Iff, const Function & x,
                const tag::LessThan &, double c1,
                const Function & v2, const tag::IfLessThan, double c2,
                const Function & v3, const tag::IfLessThan, double c3,
                const Function & v4, const tag::IfLessThan, double c4,
	              const Function & v5, const tag::Otherwise &             )
	:	arg { x }, values { v1, v2, v3, v4, v5 }, cuts { c1, c2, c3, c4 }
	{	assert ( c1 < c2 );  assert ( c1 < c2 );  assert ( c3 < c4 );  }

	inline Step ( const Function & aarg, std::vector < Function > & vals,
                std::vector < double > cts )
	: arg { aarg }, values { vals }, cuts { cts }
	{	assert ( vals.size() == cts.size() + 1 );  }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif

}; // end of  class Function::Step

//-----------------------------------------------------------------------------------------//

class Function::Vector : public Function::Core
	
{	public :

	// size_t nb_of_components ( )  stays pure virtual from Function::Core
	
	// Function component ( )  stays pure virtual from Function::Core

	virtual std::vector<double> set_value_on_cell
	( Cell::Core *, const std::vector<double> & ) = 0;
	// assign a numeric vector to the function on the cell and return that vector
	
	virtual std::vector<double> get_value_on_cell ( Cell::Core * ) const = 0;

	Function deriv ( Function ) const;
	//  virtual from Function::Core, here forbids execution, to change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core, here forbids execution
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::Aggregate : public Function::Vector
	
{	public :

	std::vector < Function > components;
	
	inline Aggregate ( const tag::ReserveSize &, size_t s )
	:	Function::Vector()
	{	components.reserve(s);  }
		
	size_t nb_of_components ( ) const;  // virtual from Function::Core, through Function::Vector

	Function component ( size_t i );
	// virtual from Function::Core, through Function::Vector
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), to change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Vector
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif
};


//-----------------------------------------------------------------------------------------//

inline Function operator&& ( Function f, Function g )

{	size_t nf = f.nb_of_components(), ng = g.nb_of_components();
	Function::Aggregate * res = new Function::Aggregate ( tag::reserve_size, nf + ng ) ;
	for ( size_t i = 0; i < nf; i++ ) res->components.emplace_back ( f[i] );
	for ( size_t i = 0; i < ng; i++ ) res->components.emplace_back ( g[i] );
	return Function ( tag::whose_core_is,	res );                                         }

//-----------------------------------------------------------------------------------------//

inline bool operator< ( const Function & f, const Function & g )
{	return f.core < g.core;  }
// needed for map 'jacobian' in class Function::Diffeomorphism

//-----------------------------------------------------------------------------------------//

class Function::Diffeomorphism : public Function::Vector
	
{	public :

	Function geom_coords, master_coords, back_geom_coords;
	// back_geom_coords are expressions involving the master coordinates
	// mathematically, they correspond to the same function geom_coords

	std::map < Function, Function > jacobian;
	// conceptually, the jacobian contains the derivatives of master_coords
	// with respect to geom_coords
	// this is implemented backwards, by computing the derivatives of
	// back_geom_coords with respect to master_coords
	// then taking the inverse matrix

	Function det;  // determinant of the jacobian matrix
	               // should be positive, so may be used directly for integration

	Diffeomorphism ( const Function & gc, const Function & mc, const Function & bmc );
	
	size_t nb_of_components ( ) const;  // virtual from Function::Core, through Function::Vector

	Function component ( size_t i );
	// virtual from Function::Core, through Function::Vector
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), to change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Vector
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function::Function ( const tag::Diffeomorphism &, const Function & geom_coords,
                            const Function & master_coords, const Function & geom_back_coords )
: Function ( tag::whose_core_is, new Function::Diffeomorphism
                ( geom_coords, master_coords, geom_back_coords ) )
{	}

//-----------------------------------------------------------------------------------------//

class Function::Composition : public Function::Scalar
	
// an expression 'base' involving master coordinates ( e.g.  1. - xi - eta )
// composed with a 'diffeom'orphism sending it in the physical space

{	public :

	Function base, diffeom;
	
	inline Composition ( const Function & b, const Function & d )
	:	Function::Scalar(), base (b), diffeom (d)
	{ }
		
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	double get_value_on_cell ( Cell::Core * ) const;  // virtual from Function::Scalar

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//

inline Function::Function
(	const Function & expr, const tag::ComposedWith &, const Function & diffeo )
// an 'expr'ession involving master coordinates ( e.g.  1. - xi - eta )
// composed with a 'diffeo'morphism sending it in the physical space
: Function ( tag::non_existent )
{	Function::Constant * exc = dynamic_cast < Function::Constant * > ( expr.core );
	if ( exc ) (*this) = expr;
	else (*this) = Function
	   ( tag::whose_core_is, new Function::Composition ( expr, diffeo ) );            }

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField
	
{	public :

  Field::Core * field;

	inline CoupledWithField ( Field::Core * f ) : field { f }
	{	assert ( f );  }
		
	class Scalar;  class Vector;
};

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Scalar
: public Function::Scalar,
	public Function::CoupledWithField
	
{	public :

	inline Scalar ( Field::Scalar * f )
	:	Function::Scalar(),
		Function::CoupledWithField ( f )
	{	assert ( f->nb_of_components() == 1 );  }
		
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	double get_value_on_cell ( Cell::Core * ) const;  // virtual from Function::Scalar

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};

//-----------------------------------------------------------------------------------------//

class Function::CoupledWithField::Vector
: public Function::Aggregate,
	public Function::CoupledWithField

// inheriting from Function::Aggregate means simply that there is a 'components' member

// quando calculamos as componentes duma funcao associada a um campo vectorial,
// nao precisamos de criar as componentes do campo vectorial
// podemos ter uma funcao escalar associada a um indice dentro dum campo vectorial
	
{	public :

	inline Vector ( Field::Block * f )
	:	Function::Aggregate ( tag::reserve_size, f->nb_of_components() ),
		Function::CoupledWithField ( f )
	{	size_t n = f->nb_of_components();
		for ( size_t j = 0; j < n; j++ )
//		components[j] = Function ( tag::whose_core_is,
			components.emplace_back ( tag::whose_core_is,
				new Function::CoupledWithField::Scalar ( field->component(j) ) );  }
		
	size_t nb_of_components ( ) const;
	// virtual from Function::Core through Function::Vector, Function::Aggregate
	
	Function component ( size_t i );
	// virtual from Function::Core, through Function::Vector, Function::Aggregate
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector, through Function::Aggregate

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), may change

	Function replace ( const Function & x, const Function & y );
	//  virtual from Function::Core, through Function::Vector
	
	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void );
	//    defined by Function::Vector (execution forbidden)
	#endif
};

//-----------------------------------------------------------------------------------------//

#ifndef NDEBUG
inline std::string Function::repr ( ) const  { return core->repr();  }
#endif

inline Function::Scalar * Function::core_to_scalar ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Scalar * res = dynamic_cast < Function::Scalar * > ( f );
	assert ( res );  return res;                                         }
#else  // no debug
	{	return static_cast < Function::Scalar * > ( f );  }
#endif
	
inline Function::Vector * Function::core_to_vector ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Vector * res = dynamic_cast < Function::Vector * > ( f );
	assert ( res );  return res;                                         }
#else  // no debug
	{	return static_cast < Function::Vector * > ( f );  }
#endif
	
inline Function::Constant * Function::core_to_constant ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Constant * res = dynamic_cast < Function::Constant * > ( f );
	assert ( res );  return res;                                         }
#else  // no debug
	{	return static_cast < Function::Constant * > ( f );  }
#endif
	
inline Function::Sum * Function::core_to_sum ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Sum * res = dynamic_cast < Function::Sum * > ( f );
	assert ( res );  return res;                                   }
#else  // no debug
	{	return static_cast < Function::Sum * > ( f );  }
#endif
	
inline Function::Product * Function::core_to_product ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Product * res = dynamic_cast < Function::Product * > ( f );
	assert ( res );  return res;                                           }
#else  // no debug
	{	return static_cast < Function::Product * > ( f );  }
#endif
	
inline Function::Diffeomorphism * Function::core_to_diffeom ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Diffeomorphism * res = dynamic_cast < Function::Diffeomorphism * > ( f );
	assert ( res );  return res;                                           }
#else  // no debug
	{	return static_cast < Function::Diffeomorphism * > ( f );  }
#endif
	
//-----------------------------------------------------------------------------------------//

Function operator+ ( const Function & f, const Function & g );
Function operator* ( const Function & f, const Function & g );
Function power ( const Function & f, double e );

inline Function operator- ( const Function & f )
{	Function minus_one ( -1. );
	return minus_one * f;         }

inline Function operator- ( const Function & f, const Function & g )
{	return f + (-g);  }

inline Function operator+= ( Function & f, const Function & g )
{	return f = f + g;  }

inline Function operator*= ( Function & f, const Function & g )
{	return f = f * g;  }

inline Function operator/ ( const Function & f, const Function & g )
{	return f * power ( g, -1. );  }

inline Function operator/= ( Function & f, const Function & g )
{	return f = f / g;  }

inline Function Function::deriv ( const Function & x ) const  // derivative with respect to x
{	return this->core->deriv ( x );  }
//-----------------------------------------------------------------------------------------//

namespace tag  {  struct Threshold { };  static const Threshold threshold;  }

inline Function sign_of ( const Function & f )
{	return Function ( tag::whose_core_is, new Function::Step
      ( Function (-1.), tag::iff, f, tag::less_than, 0., Function(1.), tag::otherwise ) );  }

inline Function smooth_abs ( const Function & f, const tag::Threshold &, double d )
// a smooth (C1) approximation of the absolute value
// 'd' is a threshold; if |f| >= d then smooth_abs(f) == |f|
// below 'd', the function is smooth but far from the absolute value (never less than d/2)
{	assert ( d > 0. );
	return Function ( tag::whose_core_is, new Function::Step
	      ( -f, tag::iff, f, tag::less_than, -d,
          ( f*f + d*d ) / (2.*d), tag::if_less_than, d, f, tag::otherwise ) );  }
//	Function abs_f = sign_of ( f ) * f;
//	return abs_f * abs_f / ( d + abs_f );  }

inline Function smooth_min
( const Function & f, const Function & g, const tag::Threshold &, double d )
// a smooth (C1) approximation of the minimum between two functions
{	return 0.5 * ( f + g - smooth_abs ( f-g, tag::threshold, d ) );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const tag::Threshold &, double d )
// a smooth (C1) approximation of the minimum between three functions
{	return smooth_min ( f, smooth_min ( g, h, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const tag::Threshold &, double d                                                )
// a smooth (C1) approximation of the minimum between four functions
{	return smooth_min ( smooth_min ( f, g, tag::threshold, d ),
                      smooth_min ( h, i, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j ,const tag::Threshold &, double d                            )
// a smooth (C1) approximation of the minimum between five functions
{	return smooth_min ( smooth_min ( f, g, h, tag::threshold, d ),
	                    smooth_min ( i, j, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_min
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j, const Function & k, const tag::Threshold &, double d        )
// a smooth (C1) approximation of the minimum between six functions
{	return smooth_min ( smooth_min ( f, g, tag::threshold, d ),
											smooth_min ( h, i, tag::threshold, d ),
	                    smooth_min ( j, k, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const tag::Threshold &, double d )
// a smooth (C1) approximation of the maximum between two functions
{	return 0.5 * ( f + g + smooth_abs ( f-g, tag::threshold, d ) );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const tag::Threshold &, double d )
// a smooth (C1) approximation of the maximum between three functions
{	return smooth_max ( f, smooth_max ( g, h, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const tag::Threshold &, double d                                                )
// a smooth (C1) approximation of the maximum between four functions
{	return smooth_max ( smooth_max ( f, g, tag::threshold, d ),
	                    smooth_max ( h, i, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j ,const tag::Threshold &, double d                            )
// a smooth (C1) approximation of the maximum between five functions
{	return smooth_max ( smooth_max ( f, g, h, tag::threshold, d ),
	                    smooth_max ( i, j, tag::threshold, d ), tag::threshold, d );  }

inline Function smooth_max
( const Function & f, const Function & g, const Function & h, const Function & i,
  const Function & j, const Function & k, const tag::Threshold &, double d        )
// a smooth (C1) approximation of the maximum between six functions
{	return smooth_max ( smooth_max ( f, g, tag::threshold, d ),
	                    smooth_max ( h, i, tag::threshold, d ),
	                    smooth_max ( j, k, tag::threshold, d ), tag::threshold, d );  }
										
	
//-----------------------------------------------------------------------------------------//

class Function::TakenOnCell	

{	public :

  Function::Core * f;
	// Function::TakenOnCell should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here
	
	Cell::Core * cll;

	TakenOnCell ( const Function::TakenOnCell & other ) = delete;
	TakenOnCell ( const Function::TakenOnCell && other ) = delete;
	
	Function::TakenOnCell & operator= ( const Function::TakenOnCell & other )
	{	if ( f->nb_of_components() == 1 )
		{	(*this) = static_cast <double> (other);
			return *this;                            }
		else
		{	(*this) = static_cast <std::vector<double>> (other);
			return *this;                                          }  }

	Function::TakenOnCell & operator= ( const Function::TakenOnCell && other )
	{	if ( f->nb_of_components() == 1 )
		{	(*this) = static_cast <double> (other);
			return *this;                            }
		else
		{	(*this) = static_cast <std::vector<double>> (other);
			return *this;                                          }  }
	
	inline operator double() const
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->get_value_on_cell ( cll );                    }
		
	inline double operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, x );                 }

	inline double operator+= ( const double & x )
	// can be used like in  f(cll) += 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) + x ); }

	inline double operator-= ( const double & x )
	// can be used like in  f(cll) -= 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) - x ); }

	inline double operator*= ( const double & x )
	// can be used like in  f(cll) *= 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) * x ); }

	inline double operator/= ( const double & x )
	// can be used like in  f(cll) /= 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, f_scalar->get_value_on_cell(cll) / x ); }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f(cll)
	{	Function::Vector * f_vect = Function::core_to_vector ( f );
	  return f_vect->get_value_on_cell ( cll );                    }
		
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Function::Vector * f_vect = Function::core_to_vector ( f );
		return f_vect->set_value_on_cell ( cll, x );                      }

};

//-----------------------------------------------------------------------------------------//

class Function::Equality

// the result of an equality comparison between Functions

{	public :

	Function lhs, rhs;
	
};

//-----------------------------------------------------------------------------------------//

inline Function::Equality operator== ( const Function & f, const Function & g )
{	return Function::Equality { f, g };  }


inline size_t Function::nb_of_components ( ) const
{	return this->core->nb_of_components();  }


inline Function Function::operator[] ( size_t i ) const
{	assert ( i < this->nb_of_components() );
	Function::Scalar * f_scalar = dynamic_cast < Function::Scalar * > ( this->core );
	if ( f_scalar )
	{	assert ( this->nb_of_components() == 1 );
		return * this;                             }
  Function::Vector * f_vector = dynamic_cast < Function::Vector * > ( this->core );
	assert ( f_vector );
	return Function ( tag::whose_core_is, f_vector->component(i).core );  }


inline Function::TakenOnCell Function::operator() ( const Cell & cll ) const
{	return Function::TakenOnCell { this->core, cll.core };   }


}  // namespace maniFEM


#endif
// ifndef MANIFEM_FUNCTION_H
