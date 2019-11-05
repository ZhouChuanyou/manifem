
// maniFEM function.cpp 2019.11.03

#include "math.h"
#include <sstream>

#include "function.h"

using namespace maniFEM;

std::shared_ptr < Function::Constant > Function::Zero ( new Function::Constant ( 0. ) );
std::shared_ptr < Function::Constant > Function::One ( new Function::Constant ( 1. ) );
std::shared_ptr < Function::Constant > Function::MinusOne ( new Function::Constant ( -1. ) );

size_t Function::Scalar::nb_of_components ( ) const  // virtual from Function::Core
{	return 1;  }

size_t Function::Aggregate::nb_of_components ( ) const
// virtual from Function::Core, through Function::Vector
{	return this->components.size();  }

size_t Function::CoupledWithField::Vector::nb_of_components ( ) const
// virtual from Function::Core through Function::Vector
{	return this->field->nb_of_components();  }


std::shared_ptr<Function::Scalar> Function::Aggregate::component ( size_t i )
// virtual from Function::Vector
{	assert ( i < this->components.size() );
	return components[i];                    }

std::shared_ptr<Function::Scalar> Function::CoupledWithField::Vector::component ( size_t i )
// virtual from Function::Vector
{	size_t n = this->field->nb_of_components();
	assert ( i < n );
	if ( this->components.size() )
	{	assert ( this->components.size() == n );
		return this->components[i];               }
	// this->components.size() == 0
	this->components.resize(n);
	for ( size_t j = 0; j < n; j++ )
		this->components[j] = std::shared_ptr < Function::Scalar >	
			( new Function::CoupledWithField::Scalar ( this->field->component(j) ) );
	return this->components[i];                                                     }
	
	
double Function::ArithmeticExpression::set_value_on_cell
( Cell::Core * cll, const double & )  // virtual from Function::Scalar
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot assign to an arithmetic expression." << std::endl;
	exit ( 1 );                                                              }
	
double Function::Constant::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ return this->val;  }

double Function::Sum::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator
		it = this->terms.begin();
	double sum = 0.;
	for ( ; it != this->terms.end(); it++ )  sum += (*it)->get_value_on_cell ( cll );
	return sum;                                                                     }
	
double Function::Product::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator
		it = this->factors.begin();
	double prod = 1.;
	for ( ; it != this->factors.end(); it++ )  prod *= (*it)->get_value_on_cell ( cll );
	return prod;                                                                      }
	
double Function::Power::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ return pow ( this->base->get_value_on_cell(cll), this->exponent );   }
	
double Function::JustTesting::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar, through Function::ArithmeticExpression
{ assert ( false );   }
	
std::vector<double> Function::Aggregate::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ size_t n = this->nb_of_components();
	assert ( n == x.size() );
	for ( size_t i = 0; i < n; i++ )
		this->components[i]->set_value_on_cell ( cll, x[i] );
	return x;                                                 }

std::vector<double> Function::Aggregate::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Vector
{ size_t n = this->nb_of_components();
	std::vector<double> result ( n );
	for ( size_t i = 0; i < n; i++ )
		result[i] = this->components[i]->get_value_on_cell ( cll );
	return result;                                                 }

double Function::CoupledWithField::Scalar::set_value_on_cell
( Cell::Core * cll, const double & x )  // virtual from Function::Scalar
{ return this->field->on_cell(cll) = x;  }

double Function::CoupledWithField::Scalar::get_value_on_cell ( Cell::Core * cll ) const
// virtual from Function::Scalar
{ return this->field->on_cell(cll).reference();  }
	
std::vector<double> Function::CoupledWithField::Vector::set_value_on_cell
( Cell::Core * cll, const std::vector<double> & x )  // virtual from Function::Vector
{ return this->field->on_cell(cll) = x;  }

std::vector<double> Function::CoupledWithField::Vector::get_value_on_cell
( Cell::Core * cll ) const  // virtual from Function::Vector
{ return this->field->on_cell(cll);  }
	

#ifndef NDEBUG	

std::string Function::Constant::repr ( const Function::From & from ) const
{	std::stringstream ss;
	ss << this->val; 
	std::string s = ss.str();
	if ( this->val < 0. ) s = "(" + s + ")";
	return s;                         }

std::string Function::Sum::repr ( const Function::From & from ) const
{ std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator
		it = this->terms.begin();
	assert ( it != terms.end() );
	std::string s = (*it)->repr ( Function::from_sum );
	for ( it++; it != terms.end(); it++ )
		s = s + '+' + (*it)->repr ( Function::from_sum );
	if ( ( from == Function::from_power ) || ( from == Function::from_product ) )
		s = "(" + s + ")";
	return s;                                                               }

std::string Function::Product::repr ( const Function::From & from ) const
{ std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator
		it = this->factors.begin();
	assert ( it != factors.end() );
	std::string s = (*it)->repr ( Function::from_product );
	for ( it++; it != factors.end(); it++ )
		s = s + '*' + (*it)->repr ( Function::from_product );
	if ( from == Function::from_power ) s = "(" + s + ")";
	return s;                                                    }

std::string Function::Power::repr ( const Function::From & from ) const
{	std::string s = base->repr ( Function::from_power ) + "^";
	std::stringstream ss;
	if ( exponent >= 0. )  ss << exponent;
	else  ss << "(" << exponent << ")";
	return s+ss.str();                      }

std::string Function::JustTesting::repr ( const Function::From & from ) const
{	return this->name;  }

std::string Function::Vector::repr ( const Function::From & from ) const
{	assert ( false );  }

std::string Function::CoupledWithField::Scalar::repr ( const Function::From & from ) const
{	return "scalar";  }

#endif // DEBUG


std::shared_ptr<Function::Core> Function::sum
( std::shared_ptr<Function::Core> f, std::shared_ptr<Function::Core> g )
// static

{	// both should be scalar :
	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f );
	std::shared_ptr < Function::Scalar > g_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( g );
	assert ( f_scalar );  assert ( g_scalar );

	// if one of them is zero :
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f );
	if ( f_const )
		if ( f_const == Function::Zero ) return g;
	std::shared_ptr < Function::Constant > g_const =
		std::dynamic_pointer_cast < Function::Constant > ( g );
	if ( g_const )
		if ( g_const == Function::Zero ) return f;
	
	// if one of them is a sum, or both :
	std::shared_ptr < Function::Sum > f_sum =
		std::dynamic_pointer_cast < Function::Sum > ( f );
	std::shared_ptr < Function::Sum > g_sum =
		std::dynamic_pointer_cast < Function::Sum > ( g );

	std::shared_ptr < Function::Sum > result { new Function::Sum };  // empty sum
	if ( g_sum )  // g is a sum
	{	std::forward_list<std::shared_ptr<Function::Scalar>>::iterator it_g;
		for ( it_g = g_sum->terms.begin(); it_g != g_sum->terms.end(); it_g++ )
			result->terms.push_front ( *it_g );                             }
	else  result->terms.push_front ( g_scalar );
	if ( f_sum )  // f is a sum
	{	std::forward_list<std::shared_ptr<Function::Scalar>>::iterator it_f;
		for ( it_f = f_sum->terms.begin(); it_f != f_sum->terms.end(); it_f++ )
			result->terms.push_front ( *it_f );                             }
	else  result->terms.push_front ( f_scalar );

	return result;                                                          }


std::shared_ptr<Function::Core> Function::difference
( std::shared_ptr<Function::Core> f, std::shared_ptr<Function::Core> g )
// static
{	return Function::sum ( f, Function::product ( Function::MinusOne, g ) );  }


std::shared_ptr<Function::Core> Function::product
( std::shared_ptr<Function::Core> f, std::shared_ptr<Function::Core> g )
// static

{	// both should be scalar
	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f );
	std::shared_ptr < Function::Scalar > g_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( g );
	assert ( f_scalar );  assert ( g_scalar );
	
	// if any one of them is zero or one :
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f );
	if ( f_const )
	{	if ( f_const == Function::Zero ) return Function::Zero;
		if ( f_const == Function::One )  return g;                }
	std::shared_ptr < Function::Constant > g_const =
		std::dynamic_pointer_cast < Function::Constant > ( g );
	if ( g_const )
	{	if ( g_const == Function::Zero ) return Function::Zero;
		if ( g_const == Function::One ) return f;                }
	
	// if one of them is a product, or both :
	std::shared_ptr < Function::Product > f_prod =
		std::dynamic_pointer_cast < Function::Product > ( f );
	std::shared_ptr < Function::Product > g_prod =
		std::dynamic_pointer_cast < Function::Product > ( g );

	std::shared_ptr < Function::Product > result { new Function::Product };  // empty product
	if ( g_prod )  // g is a product
	{	std::forward_list<std::shared_ptr<Function::Scalar>>::iterator it_g;
		for ( it_g = g_prod->factors.begin(); it_g != g_prod->factors.end(); it_g++ )
			result->factors.push_front ( *it_g );                             }
	else  result->factors.push_front ( g_scalar );
	if ( f_prod )  // f is a product
	{	std::forward_list<std::shared_ptr<Function::Scalar>>::iterator it_f;
		for ( it_f = f_prod->factors.begin(); it_f != f_prod->factors.end(); it_f++ )
			result->factors.push_front ( *it_f );                             }
	else  result->factors.push_front ( f_scalar );

	return result;                                                          }


std::shared_ptr<Function::Core> Function::division
( std::shared_ptr<Function::Core> f, std::shared_ptr<Function::Core> g )
// static
{	return Function::product ( f, Function::power ( g, -1. ) );  }


std::shared_ptr<Function::Core> Function::power ( std::shared_ptr<Function::Core> f, double e )
// static

{	// if f is a product, return a product, otherwise return a power
	// more simplifications can be done, but for the moment we don't go much deeper

	std::shared_ptr < Function::Scalar > f_scalar =
		std::dynamic_pointer_cast < Function::Scalar > ( f );
	assert ( f_scalar );
	
	std::shared_ptr < Function::Constant > f_const =
		std::dynamic_pointer_cast < Function::Constant > ( f );
	if ( f_const )
	{	if ( f_const == Function::Zero ) return Function::Zero;
		if ( f_const == Function::One )  return Function::One;
		return std::shared_ptr<Function::Core>
		{ new Function::Constant ( pow ( f_const->val, e ) ) };   }

	std::shared_ptr < Function::Product > f_prod =
		std::dynamic_pointer_cast < Function::Product > ( f );

	if ( f_prod )  // f is a product
	{	std::shared_ptr < Function::Product > result { new Function::Product };  // empty product
		std::forward_list<std::shared_ptr<Function::Scalar>>::iterator it;
		for ( it = f_prod->factors.begin(); it != f_prod->factors.end(); it++ )
		{	std::shared_ptr < Function::Scalar > g = *it;
			std::shared_ptr < Function::Scalar > g_pow =
				std::dynamic_pointer_cast < Function::Scalar > ( Function::power ( g, e ) );
			assert ( g_pow );
			result->factors.push_front ( g_pow );  }
		return result;                                                                      }

	return std::shared_ptr<Function::Core> { new Function::Power ( f_scalar, e ) };         }


std::shared_ptr<Function::Core> Function::Constant::deriv
( std::shared_ptr<Function::Core> ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	return Function::Zero;  }


std::shared_ptr<Function::Core> Function::Sum::deriv
( std::shared_ptr<Function::Core> x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	std::shared_ptr<Function::Core> result = Function::Zero;
	std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator it;
	for ( it = this->terms.begin(); it != this->terms.end(); it++ )
		result = Function::sum ( result, (*it)->deriv(x) );
	return result;                                                    }

std::shared_ptr<Function::Core> Function::Product::deriv
( std::shared_ptr<Function::Core> x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	std::shared_ptr<Function::Core> result = Function::Zero;
	std::forward_list<std::shared_ptr<Function::Scalar>>::const_iterator it1, it2;
	size_t c1, c2;
	for ( it1 = this->factors.begin(), c1 = 0;
        it1 != this->factors.end(); it1++, c1++ )
	{	std::shared_ptr<Function::Core> partial_res = Function::One;
		for ( it2 = this->factors.begin(), c2 = 0;
	        it2 != this->factors.end(); it2++, c2++ )
			if ( c1 == c2 )  // later we can eliminate c1 and c2
			{	assert ( it1 == it2 );
				partial_res = Function::product ( partial_res, (*it2)->deriv(x) );  }
			else				
			{	assert ( it1 != it2 );
				partial_res = Function::product ( partial_res, *it2 );  }
		result = Function::sum ( result, partial_res );                              }
	return result;                                                                     }

std::shared_ptr<Function::Core> Function::Power::deriv
( std::shared_ptr<Function::Core> x ) const
//  virtual from Function::Core, through Function::Scalar, Function::ArithmeticExpression
{	// we use Function constructor as a factory function for Function::Constant
	Function expo ( this->exponent );
	std::shared_ptr<Function::Core> result = expo.core;
	std::shared_ptr<Function::Core> new_power =
		Function::power ( this->base, this->exponent - 1. );
	result = Function::product ( result, new_power );
	std::shared_ptr<Function::Core> deriv_base = this->base->deriv ( x );
	result = Function::product ( result, deriv_base );
	return result;                                                           }

std::shared_ptr<Function::Core> Function::JustTesting::deriv
( std::shared_ptr<Function::Core> x ) const
{ if ( this == x.get() ) return Function::One;
	return Function::Zero;                        }

std::shared_ptr<Function::Core> Function::CoupledWithField::Scalar::deriv
( std::shared_ptr<Function::Core> x ) const
{ if ( this == x.get() ) return Function::One;
	return Function::Zero;                        }

std::shared_ptr<Function::Core> Function::Vector::deriv
( std::shared_ptr<Function::Core> x ) const
{ assert ( false );  }



