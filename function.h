
// maniFEM function.h 2019.11.23

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

class Manifold;
	

class Function

// just a a thin wrapper around Function::Core
// I guess this is what they call "delegation" ...
// useful for arithmetic operators and for vector functions

{	public :

	class Core;
	class Scalar;  class ArithmeticExpression;  class Constant;
	class Vector;  class Aggregate;  class CoupledWithField;
	class Sum;  class Product;  class Power;  class Sin;  class Cos;
	class Equality;
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

	static Function Zero, One, MinusOne;

	inline Function deriv ( const Function & x ) const;  // derivative with respect to x

	struct Catalogue
	{	std::shared_ptr < Function::Core > self;
		std::map < Function::Core *, std::shared_ptr < Function::Core > > deriv;  };
	
	#ifndef NDEBUG
	inline std::string repr ( );
	enum From { from_void, from_sum, from_product, from_power, from_function };
	#endif

	static inline Function::Scalar * core_to_scalar ( Function::Core * f );
	
	static inline Function::Vector * core_to_vector ( Function::Core * f );

};


class Function::Core

// just a base for classes like Constant, Sum, Product and many others

{	public :

	// Manifold * manifold;
	
	virtual ~Core ( ) { };

	virtual size_t nb_of_components ( ) const = 0;

	// Function component ( size_t i )
	//    not defined for Function::Scalar

	virtual Function deriv ( Function ) const = 0;

	#ifndef NDEBUG	
	virtual std::string repr ( const Function::From & from = Function::from_void ) const = 0;
	#endif
};


class Function::Scalar : public Function::Core
	
{	public :

	size_t nb_of_components ( ) const;  // virtual from Function::Core, here returns 1

	// Function component ( size_t i )  not defined
	// we should return 'this' but we need a shared_ptr
	// returning a new, fresh shared_ptr towards 'this' would be a disaster

	virtual double set_value_on_cell ( Cell::Core *, const double & ) = 0;
	// assign a numeric value to the function on the cell and return that value
	
	virtual double get_value_on_cell ( Cell::Core * ) const = 0;

	// Function deriv ( Function )
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

	// Function deriv ( Function )
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

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

	
inline Function::Function ( double c )
:	Function ( tag::whose_core_is, new Function::Constant ( c ) )
{	}


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

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


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

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


class Function::Power : public Function::ArithmeticExpression
	
{	public :

	Function base;
	double exponent;

	inline Power ( Function b, double e )
	:	base { b }, exponent { e }
	{	assert ( std::dynamic_pointer_cast < Function::Scalar > ( b.core ) );  }
		
	
	// size_t nb_of_components ( )  defined by Function::Scalar, returns 1

	// double set_value_on_cell ( Cell::Core *, const double & )
	//   defined by Function::ArithmeticExpression (execution forbidden)

	double get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Scalar, through Function::ArithmeticExpression

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


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

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};


inline Function sin ( const Function & f )
{	return Function ( tag::whose_core_is,
		new Function::Sin ( f ) );  }


	
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

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
	#endif
};

inline Function cos ( const Function & f )
{	return Function ( tag::whose_core_is,
		new Function::Cos ( f ) );  }


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

	Function deriv ( Function ) const;
	//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	// virtual from Function::Core, through Function::Scalar
	#endif
};


class Function::Vector : public Function::Core
	
{	public :

	// size_t nb_of_components ( )  stays pure virtual from Function::Core
	
	virtual Function component ( size_t i ) = 0;

	virtual std::vector<double> set_value_on_cell
	( Cell::Core *, const std::vector<double> & ) = 0;
	// assign a numeric vector to the function on the cell and return that vector
	
	virtual std::vector<double> get_value_on_cell ( Cell::Core * ) const = 0;

	Function deriv ( Function ) const;
	//  virtual from Function::Core, here forbids execution, to change

	#ifndef NDEBUG	
	std::string repr ( const Function::From & from = Function::from_void ) const;
	//  virtual from Function::Core, here forbids execution
	#endif
};


class Function::Aggregate : public Function::Vector
	
{	public :

	std::vector < Function > components;
	
	inline Aggregate ( const tag::ReserveSize &, size_t s )
	:	Function::Vector()
	{	components.reserve(s);  }
		
	size_t nb_of_components ( ) const;  // virtual from Function::Core, through Function::Vector

	Function component ( size_t i );
	// virtual from Function::Vector
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
	// virtual from Function::Vector

	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector

	// Function deriv ( Function )
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

	Function deriv ( Function ) const;
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
	// virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> set_value_on_cell ( Cell::Core *, const std::vector<double> & );
  // virtual from Function::Vector, through Function::Aggregate
	
	std::vector<double> get_value_on_cell ( Cell::Core * ) const;
	// virtual from Function::Vector, through Function::Aggregate

	// Function deriv ( Function )
	//    defined by Function::Vector (execution forbidden), may change

	#ifndef NDEBUG	
	// std::string repr ( const Function::From & from = Function::from_void );
	//    defined by Function::Vector (execution forbidden)
	#endif
};

#ifndef NDEBUG
inline std::string Function::repr ( )  { return core->repr();  }
#endif

inline Function::Scalar * Function::core_to_scalar ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Scalar * res = dynamic_cast < Function::Scalar * > ( f );
	assert ( res );
	return res;                                                           }
#else  // no debug
	{	return ( Function::Scalar * ) f;  }
#endif
	
inline Function::Vector * Function::core_to_vector ( Function::Core * f )  // static
#ifndef NDEBUG
{	Function::Vector * res = dynamic_cast < Function::Vector * > ( f );
	assert ( res );
	return res;                                                           }
#else  // no debug
	{	return ( Function::Vector * ) ss;  }
#endif
	
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


inline Function operator&& ( Function f, Function g )

{	size_t nf = f.nb_of_components(), ng = g.nb_of_components();
	Function::Aggregate * res = new Function::Aggregate ( tag::reserve_size, nf + ng ) ;
	for ( size_t i = 0; i < nf; i++ ) res->components.emplace_back ( f[i] );
	for ( size_t i = 0; i < ng; i++ ) res->components.emplace_back ( g[i] );
	return Function ( tag::whose_core_is,	res );                                         }

										
	
class Function::TakenOnCell	

{	public :

  Function::Core * f;
	// Function::TakenOnCell should only be used as temporary objects
	// they should be immediately converted to a (reference to a) double or vector<double>
	// so an ordinary pointer is OK here
	
	Cell::Core * cll;

	inline operator double() const
	// can be used like in  double x = f(cll)  or  cout << f(cll)
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->get_value_on_cell ( cll );                    }
		
	inline double operator= ( const double & x )
	// can be used like in  f(cll) = 2.0
	{	Function::Scalar * f_scalar = Function::core_to_scalar ( f );
		return f_scalar->set_value_on_cell ( cll, x );                 }

	inline operator std::vector<double>() const
	// can be used like in  vector<double> vec = f(cll)
	{	Function::Vector * f_vect = Function::core_to_vector ( f );
	  return f_vect->get_value_on_cell ( cll );                    }
		
	inline std::vector<double> operator= ( const std::vector<double> & x )
	// can be used like in  f(cll) = vec
	{	Function::Vector * f_vect = Function::core_to_vector ( f );
		return f_vect->set_value_on_cell ( cll, x );                      }

};


class Function::Equality

// the result of an equality comparison between Functions

{	public :

	Function lhs, rhs;
	
};


inline Function::Equality operator== ( const Function & f, const Function & g )
{	return Function::Equality { f, g };  }


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
	return Function ( tag::whose_core_is, f_vector->component(i).core );  }


inline Function::TakenOnCell Function::operator() ( const Cell & cll ) const
{	return Function::TakenOnCell { this->core.get(), cll.core };   }


}  // namespace maniFEM


#endif
// ifndef MANIFEM_FUNCTION_H
