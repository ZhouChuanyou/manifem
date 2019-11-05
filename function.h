
// maniFEM function.h 2019.11.04

#ifndef MANIFEM_FUNCTION_H
#define MANIFEM_FUNCTION_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <memory>
#include "assert.h"

#include "mesh.h"
#include "field.h"


namespace maniFEM {

class Manifold;
	

class Function

// just a a thin wrapper around Function::Core
// I guess this is what they call "delegation" ...
// useful for arithmetic operators and for vector functions

{	public :

	class Core;
	class Scalar;  class ArithmeticExpression;  class Constant;
	class Vector;  class Aggregate;  class CoupledWithField;
	class Sum;  class Product;  class Power;
	class Equation;
	class JustTesting;

	std::shared_ptr < Function::Core > core;

	inline Function ( const tag::NonExistent & )
	:	core ( nullptr )  { }

	inline Function ( const tag::WhoseCoreIs &, std::shared_ptr < Function::Core > c )
	:	core ( c )
	{	assert ( c );  }

	inline Function ( const tag::WhoseCoreIs &, Function::Core * c )
		:	Function ( tag::whose_core_is, std::shared_ptr < Function::Core > ( c ) )
	{  }  // use only with a fresh, new pointer

	inline Function ( const Function & f )
	:	core ( f.core )  { }

	inline Function ( const Function && f )
	:	core ( f.core )  { }  // is this the right way to do it ?

	inline Function ( const tag::HasSize &, size_t s );

	inline Function ( double c );

	inline Function & operator= ( const Function & m )
	{	core = m.core;
		return *this;   }
	
	inline Function & operator= ( const Function && m )
	{	core = m.core;  // is this the right way to do it ?
		return *this;   }
	
	inline size_t nb_of_components ( ) const;

	inline Function operator[] ( size_t );

	class TakenOnCell;
	inline Function::TakenOnCell operator() ( const Cell & cll ) const;

	static std::shared_ptr<Function::Core> sum
		( std::shared_ptr<Function::Core>, std::shared_ptr<Function::Core> );
	static std::shared_ptr<Function::Core> difference
		( std::shared_ptr<Function::Core>, std::shared_ptr<Function::Core> );
	static std::shared_ptr<Function::Core> product
		( std::shared_ptr<Function::Core>, std::shared_ptr<Function::Core> );
	static std::shared_ptr<Function::Core> division
		( std::shared_ptr<Function::Core>, std::shared_ptr<Function::Core> );
	static std::shared_ptr<Function::Core> power
		( std::shared_ptr<Function::Core>, double );

	static std::shared_ptr < Function::Constant > Zero;
	static std::shared_ptr < Function::Constant > One;
	static std::shared_ptr < Function::Constant > MinusOne;

	inline Function deriv ( const Function & x ) const;  // derivative with respect to x

	struct Catalogue
	{	std::shared_ptr < Function::Core > self;
		std::map < Function::Core *, std::shared_ptr < Function::Core > > deriv;  };
	
	#ifndef NDEBUG	
	enum From { from_void, from_sum, from_product, from_power };
	#endif
};


class Function::Core

// just a base for classes like Constant, Sum, Product and many others

{	public :

	// Manifold * manifold;
	
	virtual ~Core ( ) { };

	virtual size_t nb_of_components ( ) const = 0;

	// std::shared_ptr < Function::Scalar > component ( size_t i )
	//    not defined for Function::Scalar

	virtual std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const = 0;

	#ifndef NDEBUG	
	virtual std::string repr ( const Function::From & from = Function::from_void ) const = 0;
	#endif
};


class Function::Scalar : public Function::Core
	
{	public :

	size_t nb_of_components ( ) const;  // virtual from Function::Core, here returns 1

	// std::shared_ptr < Function::Scalar > component ( size_t i )  not defined
	// we should return 'this' but we need a shared_ptr
	// returning a new, fresh shared_ptr towards 'this' would be a disaster

	virtual double set_value_on_cell ( Cell::Core *, const double & ) = 0;
	// assign a numeric value to the function on the cell and return that value
	
	virtual double get_value_on_cell ( Cell::Core * ) const = 0;

	// std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> )
  //    stays pure virtual from Function::Core

#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
  //   stays pure virtual from Function::Core
	#endif
};


class Function::ArithmeticExpression : public Function::Scalar
	
// base for classes like Constant, Sum, Product
// any function for which 'set_value_on_cell' does not make sense

{	public :

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	double set_value_on_cell ( Cell::Core *, const double & );  // virtual from Function::Scalar
	// here execution forbidden
	
	// double get_value_on_cell ( Cell::Core * )  stays pure virtual from Function::Scalar

	// std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> )
	//   stays pure virtual from Function::Core, through Function::Scalar
	
	#ifndef NDEBUG	
	// string repr ( const Function::From & from = Function::from_void )
  //   stays pure virtual from Function::Core, through Function::Scalar
	#endif
};

	
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

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

	
class Function::Sum : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < std::shared_ptr < Function::Scalar > > terms;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


class Function::Product : public Function::ArithmeticExpression
	
{	public :

	std::forward_list < std::shared_ptr < Function::Scalar > > factors;
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


class Function::Power : public Function::ArithmeticExpression
	
{	public :

	std::shared_ptr < Function::Scalar > base;
	double exponent;

	inline Power ( std::shared_ptr < Function::Scalar > b, double e )
	:	base { b }, exponent { e }  { }
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


class Function::JustTesting : public Function::ArithmeticExpression
	
{	public :

	std::string name;

	inline JustTesting ( std::string s )
	:	name{s} { }

	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};


class Function::Vector : public Function::Core
	
{	public :

	// size_t nb_of_components ( )  stays pure virtual from Function::Core
	
	virtual std::shared_ptr < Function::Scalar > component ( size_t i ) = 0;

	virtual std::vector<double> set_value_on_cell
	( Cell::Core *, const std::vector<double> & ) = 0;
	// assign a numeric vector to the function on the cell and return that vector
	
	virtual std::vector<double> get_value_on_cell ( Cell::Core * ) const = 0;

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, here forbids execution, to change

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core, here forbids execution
	#endif
};


class Function::Aggregate : public Function::Vector
	
{	public :

	std::vector < std::shared_ptr < Function::Scalar > > components;
	
	inline Aggregate ( const tag::HasSize &, size_t s )
	:	Function::Vector(), components ( s )
	{	}
		
	size_t nb_of_components ( ) const;  // virtual from Function::Core, through Function::Vector

	std::shared_ptr < Function::Scalar > component ( size_t i );
	// virtual from Function::Vector
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> )
	//    defined by Function::Vector (execution forbidden), to change

	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void )
	//    defined by Function::Vector (execution forbidden)
	#endif
};

	
class Function::CoupledWithField
	
{	public :

  Field::Core * field;

	inline CoupledWithField ( Field::Core * f ) : field { f }
	{	assert ( f );  }
		
	class Scalar;  class Vector;
};


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

	std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> ) const;
	//  virtual from Function::Core, through Function::Scalar

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};


class Function::CoupledWithField::Vector
: public Function::Aggregate,
	public Function::CoupledWithField

// inheriting from Function::Aggregate means simply that there is a 'components' member
// initially, this vector will be empty (size zero)
// first time we compute a component, this vector will be populated with components
// later calls to 'component' will use the already-computed components
	
{	public :

	inline Vector ( Field::Block * f )
	:	Function::Aggregate ( tag::has_size, 0 ),
		Function::CoupledWithField ( f )
	{	}
		
	size_t nb_of_components ( ) const;
	// virtual from Function::Core through Function::Vector, Function::Aggregate
	
	std::shared_ptr < Function::Scalar > component ( size_t i );
	// virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
  // virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector, through Function::Aggregate

	// std::shared_ptr<Function::Core> deriv ( std::shared_ptr<Function::Core> )
	//    defined by Function::Vector (execution forbidden), may change

	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void );
	//    defined by Function::Vector (execution forbidden)
	#endif
};


inline Function::Function ( double c )
:	Function ( tag::non_existent )  // temporary empty wrapper
{	if ( c == 0. ) core = Function::Zero;
	else if ( c == 1. ) core = Function::One;
	else if ( c == -1. ) core = Function::MinusOne;
	else core.reset ( new Function::Constant ( c ) );   }


inline Function operator+ ( const Function & f, const Function & g )
{	return Function ( tag::whose_core_is, Function::sum ( f.core, g.core ) );  }

inline Function operator- ( const Function & f )
{	return Function ( tag::whose_core_is, Function::product ( f.core, Function::MinusOne ) );  }

inline Function operator- ( const Function & f, const Function & g )
{	return Function ( tag::whose_core_is, Function::difference ( f.core, g.core ) );  }

inline Function operator* ( const Function & f, const Function & g )
{	return Function ( tag::whose_core_is, Function::product ( f.core, g.core ) );  }

inline Function operator/ ( const Function & f, const Function & g )
{	return Function ( tag::whose_core_is, Function::division ( f.core, g.core ) );  }

inline Function Function::deriv ( const Function & x ) const  // derivative with respect to x
{	return Function ( tag::whose_core_is, this->core->deriv ( x.core ) );  }

	
class Function::TakenOnCell	

{	public :

  Function::Core * f;
	// Function::TakenOnCell should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here
	
	Cell::Core * cll;

	inline operator double() const
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	Function::Scalar * f_scalar = dynamic_cast<Function::Scalar*> ( f );
		assert ( f_scalar );
		return f_scalar->get_value_on_cell ( cll );                              }
		
	inline double operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	Function::Scalar * f_scalar = dynamic_cast<Function::Scalar*> ( f );
		assert ( f_scalar );
		return f_scalar->set_value_on_cell ( cll, x );                       }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f(cll)
	{	Function::Vector * f_vect = dynamic_cast<Function::Vector*> ( f );
		assert ( f_vect );
	  return f_vect->get_value_on_cell ( cll );                               }
		
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Function::Vector * f_vect = dynamic_cast<Function::Vector*> ( f );
		assert ( f_vect );
		return f_vect->set_value_on_cell ( cll, x );                      }

};


class Function::Equation

// the result of an equality comparison between Functions
// we only keep the left hand side
// the operator== will compute the difference between its two arguments

{	public :

	Function lhs;
	
};


inline Function::Equation operator== ( const Function & f, const Function & g )
{	return Function::Equation { f - g };  }


inline size_t Function::nb_of_components ( ) const
{	return this->core->nb_of_components();  }


inline Function Function::operator[] ( size_t i )
{	assert ( i < this->nb_of_components() );
	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( this->core );
	if ( f_scalar )
	{	assert ( this->nb_of_components() == 1 );
		return * this;                             }
	std::shared_ptr < Function::Vector > f_vector =
		std::dynamic_pointer_cast < Function::Vector > ( this->core );
	assert ( f_vector );
	return Function ( tag::whose_core_is, f_vector->component(i) );  }


inline Function::TakenOnCell Function::operator() ( const Cell & cll ) const
{	return Function::TakenOnCell { this->core.get(), cll.core };   }


}  // namespace maniFEM

#endif
// ifndef MANIFEM_FUNCTION_H
