
// finite-elem.h 2020.03.10

//    This file is part of maniFEM, a C++ library for meshes and finite elements on manifolds.

//    Copyright 2019, 2020 Cristian Barbarosie cristian.barbarosie@gmail.com
//    https://github.com/cristian-barbarosie/manifem

//    ManiFEM is free software: you can redistribute it and/or modify it
//    under the terms of the GNU Lesser General Public License as published
//    by the Free Software Foundation, either version 3 of the License
//    or (at your option) any later version.

//    ManiFEM is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//    See the GNU Lesser General Public License for more details.

//    You should have received a copy of the GNU Lesser General Public License
//    along with maniFEM.  If not, see <https://www.gnu.org/licenses/>.

#ifndef MANIFEM_FINITE_ELEM_H
#define MANIFEM_FINITE_ELEM_H

#include <iostream>
#include <vector>
#include <forward_list>
#include "assert.h"

#include "mesh.h"
#include "field.h"
#include "function.h"

namespace maniFEM {

namespace tag {
	struct gauss { };  static const gauss Gauss;
	enum gauss_quadrature { Q4, Q9 };
}

//-----------------------------------------------------------------------------------------//
	
class Integrator

// a wrapper for Gauss quadrature
// and perhaps other integration methods (e.g. symbolic integration)

{	public :

	class Core, Gauss;
	
	Integrator::Core * core;

	// constructor

	inline Integrator ( const tag::gauss &, const tag::gauss_quadrature & q );

	// destructor
	
	inline ~Integrator ()  {  delete core;  }
	
	// forbid copying

	Integrator ( const Integrator & ) = delete;
	Integrator ( const Integrator && ) = delete;
	Integrator & operator= ( const Integrator & ) = delete;
	Integrator & operator= ( const Integrator && ) = delete;
	
	// operator()
	
	inline double operator() ( const Function & f );
	
};  // end of  class Integrator


//-----------------------------------------------------------------------------------------//

class Integrator::Core

// a base class for Gauss quadrature
// and perhaps other integration methods (e.g. symbolic integration)

{	public :

	// constructor

	inline Core ()  { };

	// destructor

	virtual ~Core()  { };

	// operator()
	
	virtual double operator() ( const Function & f ) = 0;
	
};  // end of  class Integrator::Core

//-----------------------------------------------------------------------------------------//

class Integrator::Gauss : public Integrator::Core

// a base class for Gauss quadrature

{	public :

	// constructor

	inline Gauss () : Integrator::Core  { };

	// operator()
	
	double operator() ( const Function & f );  // virtual from Integrator::Core
	
};  // end of  class Integrator::Gauss

//-----------------------------------------------------------------------------------------//



inline Integrator::Integrator ( const tag::gauss &, const tag::gauss_quadrature & q )
:	core { new Integrator::Gauss ( q ) }  { }
	
inline double Integrator::operator() ( const Function & f )
{	return (*core) ( f );  }


//-----------------------------------------------------------------------------------------//

// Again ...
// What is a finite element ?


class FiniteElement

// wrapper for different types of finite elements	
	
{	public :

	class Core, withMaster, Lagrange;
	
	FiniteElement::Core * core;

	// constructor

	inline FiniteElement ( const tag::gauss &, const tag::gauss_quadrature & q );

	// destructor
	
	inline ~FiniteElement ()  {  delete core;  }
	
	// forbid copying

	FiniteElement ( const FiniteElement & ) = delete;
	FiniteElement ( const FiniteElement && ) = delete;
	FiniteElement & operator= ( const FiniteElement & ) = delete;
	FiniteElement & operator= ( const FiniteElement && ) = delete;
	
};  // end of  class FiniteElement
	
//-----------------------------------------------------------------------------------------//

class FiniteElement::Core

// a base class for FiniteElements

{	public :

	// constructor

	inline Core ()  { };

	// destructor

	virtual ~Core()  { };

	
};  // end of  class FiniteElement::Core

//-----------------------------------------------------------------------------------------//

class FiniteElement::withMaster : public FiniteElement::Core

// finite elements which use a master element

{	public :

	std::list < Function > master_coords, composed_space_coords;

	FunctionOnMesh::Mapping * transf;

	// constructor

	inline withMaster () : FiniteElement::Core ()  { };
	
};  // end of  class FiniteElement::withMaster

//-----------------------------------------------------------------------------------------//


}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_FINITE_ELEM_H
