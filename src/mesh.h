
// mesh.h 2021.05.01

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


#ifndef MANIFEM_MESH_H
#define MANIFEM_MESH_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include "assert.h"

namespace maniFEM {

// tags are used to distinguish between different versions of the same function
// or method or constructor
// they create a sort of spoken language ...
	
namespace tag {  // see paragraph 9.2 in the manual

	// tags are listed in no particular order
	
	struct IsNegative { };  static const IsNegative is_negative;
	struct IsPositive { };  static const IsPositive is_positive;
	struct ThisMeshIsPositive { };  static const ThisMeshIsPositive this_mesh_is_positive;
	struct ThisMeshIsNegative { };  static const ThisMeshIsNegative this_mesh_is_negative;
	struct RequireOrder { };  static const RequireOrder require_order;
	struct ReverseOrder { };  static const ReverseOrder reverse_order;
	struct ReverseOrderIfAny { };  static const ReverseOrderIfAny reverse_order_if_any;
	struct ReverseOf { };  static const ReverseOf reverse_of;
	struct Reversed { };  static const Reversed reversed;
	struct ReverseEachCell { };  static const ReverseEachCell reverse_each_cell;
	struct AsTheyAre { };  static const AsTheyAre as_they_are;
	struct NonExistent { };  static const NonExistent non_existent;
	struct BuildIfNotExists { };  static const BuildIfNotExists build_if_not_exists;
	struct SeenFrom { };  static const SeenFrom seen_from;
	struct Fuzzy { };  static const Fuzzy fuzzy;
	struct MayNotExist { };  static const MayNotExist may_not_exist;
	struct DoNotBother { };  static const DoNotBother do_not_bother;
	struct SurelyExists { };  static const SurelyExists surely_exists;
	struct CellsSurelyExist { };  static const CellsSurelyExist cells_surely_exist;
	struct OfDimension { };  static const OfDimension of_dim;
	                         static const OfDimension of_dimension;
	struct OfMaxDim { };  static const OfMaxDim of_max_dim;
	struct MinusOne { };  static const MinusOne minus_one;
	struct GreaterThanOne { };  static const GreaterThanOne greater_than_one;
	struct MightBeOne { };  static const MightBeOne might_be_one;
	struct Oriented { };  static const Oriented oriented;
	struct NotOriented { };  static const NotOriented not_oriented;
	struct DeepCopyOf { };  static const DeepCopyOf deep_copy_of;
	struct BuildCellsIfNec { };  static const BuildCellsIfNec build_cells_if_necessary;
	struct Progressive { };  static const Progressive progressive;
	struct StartAt { };  static const StartAt start_at;
	struct StopAt { };  static const StopAt stop_at;
	struct Towards { };  static const Towards towards;
	struct Boundary { };  static const Boundary boundary;
	struct BoundaryOf { };  static const BoundaryOf boundary_of;
	struct Bizarre { };  static const Bizarre bizarre;
	struct SizeMeshes { };  static const SizeMeshes size_meshes;
	struct BehindFace { };  static const BehindFace behind_face;
	struct InFrontOfFace { };  static const InFrontOfFace in_front_of_face;
	struct WithinMesh { };  static const WithinMesh within_mesh;
	struct Vertices { };  static const Vertices vertices;
	struct Segments { };  static const Segments segments;
	struct OverVertices { };  static const OverVertices over_vertices;
	struct OverSegments { };  static const OverSegments over_segments;
	struct CellsOfDim { };  static const CellsOfDim cells_of_dim;
	struct CellsOfMaxDim { };  static const CellsOfMaxDim cells_of_max_dim;
	struct OverCells { };  static const OverCells over_cells;
	struct OverCellsOfDim { };  static const OverCellsOfDim over_cells_of_dim;
	struct OverCellsOfMaxDim { };  static const OverCellsOfMaxDim over_cells_of_max_dim;
	struct OverCellsOfReverseOf { }; static const OverCellsOfReverseOf over_cells_of_reverse_of;
	struct Vertex { };  static const Vertex vertex; static const Vertex point;
	struct Segment { };  static const Segment segment;
	struct DividedIn { };  static const DividedIn divided_in;
	struct Triangle { };  static const Triangle triangle;
	struct Quadrangle { };  static const Quadrangle quadrangle;
	                        static const Quadrangle rectangle;
	                        static const Quadrangle square;
	                        static const Quadrangle quadrilateral;
	struct WhoseBoundaryIs { };  static const WhoseBoundaryIs whose_bdry_is;
	                             static const WhoseBoundaryIs whose_boundary_is;
	struct WhoseCoreIs { };  static const WhoseCoreIs whose_core_is;
	struct ForcePositive { };  static const ForcePositive force_positive;
	struct HasSize { };  static const HasSize has_size;
	struct ReserveSize { };  static const ReserveSize reserve_size;
	struct Pretty { };  static const Pretty pretty;
	struct Adapt { };  static const Adapt adapt;
	struct OfDegree { };  static const OfDegree of_degree;
	struct MeshIsBdry { };  static const MeshIsBdry mesh_is_bdry;
	struct MeshIsNotBdry { };  static const MeshIsNotBdry mesh_is_not_bdry;
	enum WithTriangles { with_triangles, not_with_triangles };
	struct Join { };  static const Join join;
	enum KeyForHook { tangent_vector, normal_vector, node_in_cloud };
	struct Onto { };  static const Onto onto;
	struct EntireManifold { };  static const EntireManifold entire_manifold;
	struct DesiredLength { };  static const DesiredLength desired_length;
	struct IntrinsicOrientation { };  static const IntrinsicOrientation intrinsic_orientation;
	struct InherentOrientation { };  static const InherentOrientation inherent_orientation;
	struct RandomOrientation { };  static const RandomOrientation random_orientation;
	struct Util
	{ template < class T > class Wrapper;
		class Core;
		class CellCore;  class MeshCore;
		inline static size_t assert_diff ( const size_t a, const size_t b )
		{	assert ( a >= b );  return  a - b;  }
		template < typename X, typename Y > inline static Y assert_cast ( X x )
		#ifndef NDEBUG
		{	Y y = dynamic_cast < Y > (x);  assert (y);  return y;  }
		#else
		{	Y y = static_cast < Y > (x);  return y;  }
		#endif
	};
	struct MayBeNull { };  static const MayBeNull may_be_null;
	struct FreshlyCreated { };  static const FreshlyCreated freshly_created;
	struct PreviouslyExisting { };  static const PreviouslyExisting previously_existing;
	struct ZeroWrappers { };  static const ZeroWrappers zero_wrappers;
	struct OldCore { };  static const OldCore old_core;
	struct SurelyNotNull { };  static const SurelyNotNull surely_not_null;
	struct NoNeedToDisposeOf { };  static const NoNeedToDisposeOf no_need_to_dispose_of;
	struct Empty { };  static const Empty empty;
	struct OneDummyWrapper { };  static const OneDummyWrapper one_dummy_wrapper;

}  // end of namespace tag

//-----------------------------------------------------------------------------//


// wrappers act like shared pointers towards cores
// you may want to forbid copying for your cores

template < class T >	class tag::Util::Wrapper 

{	public :

	T * core;

	inline Wrapper ( const tag::Empty & ) : core { nullptr }  { }

	inline Wrapper ( T* c, const tag::FreshlyCreated & ) : core { c }
	// c just constructed with nb_of_wrappers == 1, no need to increment
	{	assert ( c );  }
	
	inline Wrapper ( T* c, const tag::PreviouslyExisting & ) : core { c }
	{	assert ( c );
		c->nb_of_wrappers++;  }
	
	inline ~Wrapper ( )
	{	dispose_core ( tag::may_be_null );  }

	inline Wrapper operator= ( const Wrapper & w )
	{	dispose_core ( tag::may_be_null );
		set_core ( w.core, tag::may_be_null );  }
	
	inline void dispose_core ( const tag::MayBeNull & )
	{	if ( this->core )
			if ( this->core->dispose() )
				delete this->core;          }
	
	inline void dispose_core ( const tag::SurelyExists & )
	{	assert ( this->core );
		if ( this->core->dispose() )
			delete this->core;          }
	
	inline void set_core ( T* c, const tag::FreshlyCreated & )
	// c just constructed with nb_of_wrappers == 1, no need to increment
	{	dispose_core ( tag::may_be_null );
		// old core may be null
		assert ( c );
		core = c;                               }
	
	inline void set_core ( T* c, const tag::PreviouslyExisting & )
	{	dispose_core ( tag::may_be_null );
		// old core may be null
		assert ( c );
		core = c;
		c->nb_of_wrappers++;  }
	
	inline void set_core ( T* c, const tag::MayBeNull & )
	{	dispose_core ( tag::may_be_null );
		// old core may be null
		core = c;
		if ( c ) c->nb_of_wrappers++;  }
	
	inline void set_core
	( T* c, const tag::FreshlyCreated &,
	  const tag::OldCore &, const tag::SurelyNotNull & )
	// c just constructed with nb_of_wrappers == 1, no need to increment
	{	dispose_core ( tag::surely_not_null );
		assert ( c );
		core = c;                               }
	
	inline void set_core
	( T* c, const tag::PreviouslyExisting &,
	  const tag::OldCore &, const tag::SurelyNotNull & )
	{	dispose_core ( tag::surely_not_null );
		assert ( this->core );
		core = c;
		c->nb_of_wrappers++;  }
	
	inline void set_core
	( T* c, const tag::FreshlyCreated &,
	  const tag::NoNeedToDisposeOf &, const tag::OldCore & )
	// c just constructed with nb_of_wrappers == 1, no need to increment
	{	core = c; }
	
	inline void set_core
	( T* c, const tag::PreviouslyExisting &,
	  const tag::NoNeedToDisposeOf &, const tag::OldCore & )
	{	core = c;
		c->nb_of_wrappers++;  }

	class Inactive;
};


// class tag::Util::Core is intended to be used by constructing
// your own Core class which inherits from tag::Util::Core
// you will probably implement polymorphic cores
// if you only declare pointers to your own core classes, you can keep the destructor inline
// however, if you want to do things like
//   tag::Util::Core * c = pointer_to_your_own_type_of_core;
//   ... some code ...
//   delete c;
// then you must declare the destructor virtual
// see e.g. item 7 in the book of Scott Meyers, Effective C++

// cores are built with nb_of_wrappers = 1
// although they did not meet any wrapper yet (they are newborn)
// we assume that the 'new' command was issued by a wrapper
// who is waiting for them to see the light of day and is eager to embrace them

class tag::Util::Core

{	public :

	size_t nb_of_wrappers;

	inline Core ( const tag::OneDummyWrapper & ) : nb_of_wrappers { 1 } { }
	// a new core is built through a wrapper
	// thus, we initialize with one instead of zero

	inline Core ( const tag::ZeroWrappers & ) : nb_of_wrappers { 0 } { }
	// in rare cases, a core is built without a wrapper
	// like for negative cells
	
	inline ~Core ( ) { }

  inline bool dispose ( )
	{	assert ( nb_of_wrappers > 0 );
		nb_of_wrappers--;
		return nb_of_wrappers == 0;    }

	class DelegateDispose;
	class Inactive;
	static bool default_dispose ( tag::Util::Core::DelegateDispose * );
	static bool dispose_cell_with_reverse ( tag::Util::Core::DelegateDispose * );

};


class tag::Util::Core::DelegateDispose : public tag::Util::Core

// some classes need specialized 'dispose' method, e.g. cells having reverse

{	public :

	bool (*dispose_p) ( tag::Util::Core::DelegateDispose * ) { & tag::Util::Core::default_dispose };

	inline DelegateDispose ( const tag::OneDummyWrapper & )
	// a new core is built through a wrapper; thus, we initialize with one instead of zero
	:	tag::Util::Core ( tag::one_dummy_wrapper ) { };

	inline DelegateDispose ( const tag::ZeroWrappers & )
	// in rare cases, a core is built without a wrapper, like for negative cells
	:	tag::Util::Core ( tag::zero_wrappers ) { };
	
	inline ~DelegateDispose ( ) { };

	inline bool dispose ( )
	{	return dispose_p ( this );  }

};


template < class T >	class tag::Util::Wrapper<T>::Inactive
// does nothing
{	public :
	T * core;
	inline Inactive ( const tag::Empty & ) : core { nullptr } { }
	inline Inactive ( T* c, const tag::FreshlyCreated & ) : core { c } { }
	inline Inactive ( T* c, const tag::PreviouslyExisting & ) : core { c } { }
	inline ~Inactive ( ) { }
	inline void set_core
	( T* c, const tag::FreshlyCreated &,
	  const tag::NoNeedToDisposeOf &, const tag::OldCore & )  { core = c; }
};

class tag::Util::Core::Inactive
// does nothing
{	public :
	inline Inactive ( const tag::OneDummyWrapper & ) { }
	inline Inactive ( const tag::ZeroWrappers & ) { }
	inline ~Inactive ( ) { }
};

//-----------------------------------------------------------------------------//

class Cell;  class Mesh;
class Manifold;  class Function;

//-----------------------------------------------------------------------------//


class CellIterator

// a thin wrapper around a CellIterator::Core, with most methods delegated to 'core'

// iterates over all cells of a given mesh (cells of given dimension)

// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative (tag::force_positive overrides this)
// for lower dimension, returns positive cells

// or, iterates over all cells above a given cell

// see paragraphs 8.5 and 8.6 in the manual
	
{	public :

	class Core;
	
	std::unique_ptr < CellIterator::Core > core;

	inline CellIterator ( const CellIterator & it );
	inline CellIterator ( CellIterator && it );

	inline CellIterator & operator= ( const CellIterator & it );
	inline CellIterator & operator= ( const CellIterator && it );

	inline CellIterator ( const tag::WhoseCoreIs &, CellIterator::Core * c )
	:	core { c }
	{	assert ( c );  }
	
	inline ~CellIterator ( ) { };
	
	inline void reset ( );	
	inline void reset ( Cell & );	
	inline Cell operator* ( );
	inline CellIterator & operator++ ( );
	inline CellIterator & operator++ ( int );
	inline CellIterator & advance ( );
	inline bool in_range ( );

	struct Over  {  class TwoVerticesOfSeg;  class CellsOfFuzzyMesh;  };
	struct Adaptor  {  class ForcePositive;  };

};  // end of class CellIterator


//-----------------------------------------------------------------------------//
//-----------------  wrappers Cell and Mesh  ----------------------------------//
//-----------------------------------------------------------------------------//


// a cell of dimension zero is a point, see class Cell::PositiveVertex and NegativeVertex
// a cell of dimension one is a segment, see class Cell::Positive::Segment and NegativeSegment
// for cells of dimension two or more, see class Cell::Positive::HighDim and Cell::Negative::HighDim
// a cell of dimension two may be a triangle, a quadrangle or some other polygon
// a cell of dimension three may be a tetrahedron, a cube or some other polyhedron
// cells of dimension four or higher may be constructed,
// but their usefulness is questionable

// cells may be positively or negatively oriented
// see class Cell::Positive and class Cell::Negative
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// a cell is mainly defined by its boundary (which is a mesh of lower dimension)
// the orientation of a cell is nothing more than an orientation of its boundary
// see the comments on orientation in class Mesh below

// see paragraphs 1.2, 8.1, 8.7 and 9.3 in the manual

// I would very much prefer the name 'Cell::Core' instead of 'tag::Util::CellCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class Cell : public tag::Util::Wrapper < tag::Util::CellCore >

#else  // no MANIFEM_COLLECT_CM

class Cell : public tag::Util::Wrapper < tag::Util::CellCore > ::Inactive 

#endif  // MANIFEM_COLLECT_CM	

// a thin wrapper around a tag::Util::CellCore, with most methods delegated to 'core'

{	public :

	typedef tag::Util::CellCore Core;

	// Cell::Core * core  inherited from tag::Util::Wrapper < Cell::Core >

	inline Cell ( const tag::NonExistent & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( tag::empty )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( tag::empty )
	#endif  // MANIFEM_COLLECT_CM	
	{ }

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c, const tag::FreshlyCreated & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( c, tag::freshly_created )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( c, tag::freshly_created )
	#endif  // MANIFEM_COLLECT_CM	
	{	}

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c, const tag::PreviouslyExisting & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ( c, tag::previously_existing )
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Cell::Core > ::Inactive ( c, tag::previously_existing )
	#endif  // MANIFEM_COLLECT_CM	
	{	}

	inline Cell ( const tag::WhoseBoundaryIs &, Mesh & );
	inline Cell ( const tag::Vertex &, const tag::IsPositive & ispos = tag::is_positive );
	inline Cell ( const tag::Segment &, const Cell & A, const Cell & B );
	inline Cell ( const tag::Bizarre &, const Cell &, const Cell & );
	inline Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA );
	inline Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                         const Cell & CD, const Cell & DA );

	// redefine operator= !!
	
	// we are still in class Cell
	
	// methods delegated to 'core'

	inline Cell reverse
	( const tag::BuildIfNotExists & build = tag::build_if_not_exists ) const;
	inline Cell reverse ( const tag::SurelyExists & ) const;
	inline Cell reverse ( const tag::MayNotExist & ) const;
	inline Mesh boundary () const;
	inline bool exists () const  { return core != nullptr;  }
	inline bool is_positive () const;
	inline Cell get_positive ();
	inline size_t dim () const;
	inline bool has_reverse () const;
	inline Cell tip () const;
	inline Cell base () const;

	inline bool belongs_to ( const Mesh & msh, const tag::Oriented & ) const;
	inline bool belongs_to ( const Mesh & msh, const tag::NotOriented & ) const;

	// method 'glue_on_bdry_of' is intensively used when building a mesh
	// it glues 'this' cell to the boundary of 'cll'
	inline void glue_on_bdry_of ( Cell & cll );
	
	// method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll' - used mainly in remeshing
	inline void cut_from_bdry_of ( Cell & cll );
	
	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead
	inline void add_to ( Mesh & msh );
	inline void remove_from ( Mesh & msh );

	// we are still in class Cell

	inline void project ( ) const;
	
	inline void project ( const tag::Onto &, const Manifold m ) const;
		
#ifndef NDEBUG
	inline void print_everything ( );
#endif

	void print_coords();

	static std::vector < size_t > double_heap_size_pos, double_heap_size_neg,
		size_t_heap_size_pos, size_t_heap_size_neg, short_int_heap_size_pos, short_int_heap_size_neg;

	struct field_to_meshes
	{	short int counter_pos;
		short int counter_neg;
		std::list<Cell>::iterator where;
		inline field_to_meshes ( short int i, short int j )
		:	counter_pos {i}, counter_neg {j} { }
		inline field_to_meshes ( short int i, short int j,
		                         std::list<Cell>::iterator w )
		:	counter_pos {i}, counter_neg {j}, where {w} { }              };

	struct field_to_meshes_same_dim
	{	short int sign;
		std::list<Cell>::iterator where;
		inline field_to_meshes_same_dim ( short int i ) : sign {i} { }
		inline field_to_meshes_same_dim
		( short int i, std::list<Cell>::iterator w )
		:	sign {i}, where {w} { }                                       };

  class Positive;  class Negative;
	
	// any way to rewrite the names below as Positive::Vertex and such ?
	// I get an error because Cell::Core::add_to_seg expects a PositiveSegment ...
	class PositiveVertex;  class PositiveNotVertex;
	class PositiveSegment;  class PositiveHighDim;
	class NegativeVertex;  class NegativeNotVertex;
	class NegativeSegment;  class NegativeHighDim;

}; // end of  class Cell


inline bool operator== ( const Cell & c1, const Cell & c2 )
{	return c1.core == c2.core;  }

inline bool operator!= ( const Cell & c1, const Cell & c2 )
{	return c1.core != c2.core;  }

inline bool operator< ( const Cell & c1, const Cell & c2 )
{	return c1.core < c2.core;  }

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


// roughly speaking, a mesh is a collection of cells of the same dimension

// an orientation of the mesh
// is nothing more than an orientation of each of its cells (of maximum dimension)
// but these orientations cannot be arbitrary, they must be compatible
// in the sense that a face common to two cells must be seen as positive
// from one of the cells and as negative from the other cell

// as a consequence, for connected meshes there are only two possible orientations
// although nothing prevents a mesh to be disconnected, we only allow for two orientations
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// negative meshes will appear mostly as boundaries of negative cells
// that's why we do not store their core in the computer's memory
// wrappers for negative meshes are built on-the-fly, e.g. in method Cell::boundary
// their core points to the (positive) reverse, the only difference is in 'sign'

// see paragraphs 1.2, 8.1, 8.7 and 9.3 in the manual

// I would very much prefer the name 'Mesh::Core' instead of 'tag::Util::MeshCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class Mesh : public tag::Util::Wrapper < tag::Util::MeshCore >

#else  // no MANIFEM_COLLECT_CM	

class Mesh : public tag::Util::Wrapper < tag::Util::MeshCore > ::Inactive

#endif  // MANIFEM_COLLECT_CM	

// a thin wrapper around a Mesh::Core, with most methods delegated to 'core'

{	public :
	
	typedef tag::Util::MeshCore Core;

	// Mesh::Core * core  inherited from tag::Util::Wrapper < Mesh::Core >
	// there are no cores for negative meshes
	// instead, we keep here a pointer to the direct (positive) core
	// thus, for meshes, 'core' points always to a positive Mesh::Core

	// this is how we distinguish between a positive mesh and a negative one
	bool (*is_pos) ();

	// we keep here the topological dimension of the largest mesh we intend to build
	// see method 'set_max_dim' and paragraph 9.5 in the manual
	static size_t maximum_dimension_plus_one;

	// low-level constructors :

	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::FreshlyCreated &,
	              const tag::IsPositive & ispos = tag::is_positive                       )
	// c is freshly created with tag::one_dummy-wrapper
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( c, tag::freshly_created ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c, tag::freshly_created ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & Mesh::return_true }
	{	}
	
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                const tag::IsPositive & ispos = tag::is_positive                           )
	// c is either pre-existent or is freshly created with tag::zero_wrappers
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( c, tag::previously_existing ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c, tag::previously_existing ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & Mesh::return_true }
	{	}
	
	// build a negative mesh from a positive one
	// without worrying whether reverse cells exist or not
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::DoNotBother &                          )
	#ifdef MANIFEM_COLLECT_CM	
	: tag::Util::Wrapper < Mesh::Core > ( c, tag::previously_existing ),
	#else  // no MANIFEM_COLLECT_CM	
	:	tag::Util::Wrapper < Mesh::Core > ::Inactive ( c, tag::previously_existing ),
	#endif  // MANIFEM_COLLECT_CM	
		is_pos { & Mesh::return_false }
	{	}
	
	// more elaborate (high-level) constructors :
	// they call one of the above, then manipulate the mesh
	
	// build a negative mesh from a positive one, assuming all cells have reverse :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::CellsSurelyExist &                   );
	
	// build a negative mesh from a positive one, creating reverse cells if necessary :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::PreviouslyExisting &,
	              const tag::IsNegative &, const tag::BuildCellsIfNec &                    );

	inline Mesh ( const tag::DeepCopyOf &, const Mesh & msh );
	inline Mesh ( const tag::DeepCopyOf &, const Mesh & msh, const tag::Fuzzy & );

	inline Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t dim,
                const tag::IsPositive & ispos = tag::is_positive                );

	inline Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t dim,
                const tag::MinusOne &, const tag::IsPositive & ispos = tag::is_positive );

	// we are still in class Mesh
	
	inline Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
                const tag::DividedIn &, const size_t n                      );
	// builds a chain of n segment cells
	
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	// builds a triangular mesh from three sides
	// the number of divisions defined by the divisions of the sides (must be the same)
	
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         );
	// builds a rectangular mesh from four sides
	// the number of divisions are already in the sides (must be the same for opposite sides)
	// if last argument is true, each rectangular cell will be cut in two triangles

	inline Mesh
	( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	  const Cell & NE, const Cell & NW, const size_t m, const size_t n,
	  const tag::WithTriangles & wt = tag::not_with_triangles           );
	// builds a rectangular mesh from four vertices
	// if last argument is true, each rectangular cell will be cut in two triangles

	inline Mesh ( const tag::Join &, const std::list<Mesh> & l );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &,
	              const Mesh &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );

	template < typename container >
	inline Mesh ( const tag::Join &, const container & l );

	// we are still in class Mesh
	// constructors with tag::Progressive are defined in progressive.cpp
	
	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
         const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length,
	       const tag::RandomOrientation &                                                 );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
	       const tag::DesiredLength &, const Function & length, const tag::RandomOrientation & );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length,
	       const tag::InherentOrientation &                                               );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length              );

	inline Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	              const tag::DesiredLength &, const Function & length,
	              const tag::RandomOrientation &                                   )
	:	Mesh ( tag::progressive, tag::boundary, interface, tag::desired_length, length )  { }

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length, const tag::IntrinsicOrientation & );
 
	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> normal,
	       const tag::DesiredLength &, const Function & length               );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> tangent,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                 );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> tangent,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length, const tag::RandomOrientation & );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::StopAt &, const Cell & stop,
	              const tag::DesiredLength &, const Function & length,
	              const tag::RandomOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, stop,
		       tag::desired_length, length                                 )  { }

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::DesiredLength &, const Function & length,
	              const tag::InherentOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, start,
		       tag::desired_length, length, tag::inherent_orientation       )  { }

	// we are still in class Mesh
	
	void pretty_constructor ( const tag::Segment &, const Cell & A, const Cell & B,
	                          const tag::DividedIn &, const size_t n );
	void pretty_constructor ( const tag::Triangle &, const Mesh & AB,
	                          const Mesh & BC, const Mesh & CA         );
	void pretty_constructor ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, const tag::WithTriangles & wt        );

	inline void copy_all_cells_to ( Mesh & msh ) const;

	inline bool is_positive () const  {  return is_pos();  }
	inline size_t dim () const;
	inline Mesh reverse () const;

	inline size_t number_of ( const tag::Vertices & ) const;
	inline size_t number_of ( const tag::Segments & ) const;
	inline size_t number_of ( const tag::CellsOfDim &, const size_t d ) const;
	inline size_t number_of ( const tag::CellsOfMaxDim & ) const;

	inline Cell first_vertex ( ) const;
	inline Cell last_vertex ( ) const;
	inline Cell first_segment ( ) const;
	inline Cell last_segment ( ) const;

	// we are still in class Mesh
	
	inline Cell cell_in_front_of
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell cell_behind ( const Cell face, const tag::MayNotExist & ) const;

	// we are still in class Mesh
	// the eight methods below are only relevant for STSI meshes

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighb,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighb,
	  const tag::MayNotExist &                                                       ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighb,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighb,
	  const tag::MayNotExist &                                                       ) const;

	// we are still in class Mesh
	
	template < typename container >
	void join_list ( const container & l );

	// method baricenter defined in manifold.h
	void inline baricenter ( const Cell & ver, const Cell & seg );

	// iterators defined in iterator.h
	inline CellIterator iterator ( const tag::OverVertices & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrderIfAny &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &      ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &      ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &    ) const;

	// we are still in class Mesh
	
	inline CellIterator iterator ( const tag::OverSegments & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrderIfAny &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &      ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &      ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &    ) const;

	// we are still in class Mesh
	
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	 const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &                                             ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder &,
	  const tag::ForcePositive &                                           ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &                                             ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrder &,
	  const tag::ForcePositive &                                           ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &                                        ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrderIfAny &,
	  const tag::ForcePositive &                                                ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t,
		const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &                     ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &               ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &                     ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &               ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &                ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &                    ) const;

	// we are still in class Mesh
	
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &                                         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::RequireOrder &                                      ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &, const tag::RequireOrder &              ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::RequireOrder &, const tag::AsTheyAre &              ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrder &                                      ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &, const tag::ReverseOrder &              ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrder &, const tag::AsTheyAre &              ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrderIfAny &                                 ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::AsTheyAre &, const tag::ReverseOrderIfAny &         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrderIfAny &, const tag::AsTheyAre &         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ForcePositive &                                     ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ForcePositive &, const tag::RequireOrder &          ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::RequireOrder &, const tag::ForcePositive &          ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ForcePositive &, const tag::ReverseOrder &          ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrder &, const tag::ForcePositive &          ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ForcePositive &, const tag::ReverseOrderIfAny &     ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrderIfAny &, const tag::ForcePositive &     ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
		const tag::ReverseEachCell &, const tag::DoNotBother &         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::RequireOrder &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrder &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfDimension &, const size_t,
	  const tag::ReverseOrderIfAny &, const tag::ReverseEachCell &, const tag::DoNotBother & ) const;

	// we are still in class Mesh
	
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator ( const tag::OverCellsOfMaxDim &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &                                   ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
	  const tag::ForcePositive &                                 ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &                                   ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrder &,
	  const tag::ForcePositive &                                 ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &                              ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ForcePositive &                                      ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &,
		const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &           ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &     ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &           ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &     ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &      ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &          ) const;

	// we are still in class Mesh
	
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::RequireOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrder &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::AsTheyAre &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrderIfAny &, const tag::AsTheyAre & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::RequireOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::RequireOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ForcePositive &, const tag::ReverseOrderIfAny & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
	  const tag::ReverseOrderIfAny &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &,
		const tag::ReverseEachCell &, const tag::DoNotBother & ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &                         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::RequireOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &                   ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &                         ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrder &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &                   ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &                    ) const;
	inline CellIterator iterator
	( const tag::OverCells &, const tag::OfMaxDim &, const tag::ReverseOrderIfAny &,
	  const tag::ReverseEachCell &, const tag::DoNotBother &                        ) const;

	// we are still in class Mesh
	// methods draw_ps and export_msh defined in global.cpp

	void draw_ps ( std::string file_name );
	void draw_ps_3d ( std::string file_name );
	void export_msh ( std::string f, std::map<Cell::Core*,size_t> & ver_numbering );
	void export_msh ( std::string f );
	
#ifndef NDEBUG
	inline void print_everything ( );
#endif
	
	inline static void set_max_dim ( const size_t d )
	// see paragraph 9.5 in the manual
	{	maximum_dimension_plus_one = d + 1;
		Cell::double_heap_size_pos.resize ( maximum_dimension_plus_one, 0. );
		Cell::double_heap_size_neg.resize ( maximum_dimension_plus_one, 0. );
		Cell::size_t_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::size_t_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );  }

	struct Connected  {  class OneDim;  class HighDim;  };
	struct MultiplyConnected  {  class OneDim;  class HighDim; };
	class ZeroDim;  class Fuzzy;  class STSI;
	
	static bool return_true ( );
	static bool return_false ( );
	
}; // end of  class Mesh


inline bool operator== ( const Mesh & m1, const Mesh & m2 )
{	return m1.core == m2.core;  }


//-----------------------------------------------------------------------------//
//----------------      core Cells     ----------------------------------------//
//-----------------------------------------------------------------------------//


// I would very much prefer the name 'class Cell::Core' instead of 'tag::Util::CellCore'
// but it was not possible ...


#ifdef MANIFEM_COLLECT_CM	

class tag::Util::CellCore : public tag::Util::Core::DelegateDispose

#else  // no MANIFEM_COLLECT_CM	

class tag::Util::CellCore : public tag::Util::Core::Inactive

#endif  // MANIFEM_COLLECT_CM	

// abstract class
// specialized in Cell::Positive::{Vertex,Segment,HighDim} and
//                Cell::Negative::{Vertex,Segment,HighDim}

{	public :
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM
	
	// we keep numeric values here :
	std::vector < double > double_heap;
	std::vector < size_t > size_t_heap;
	std::vector < short int > short_int_heap;
	// and heterogeneous information here :
	std::map < tag::KeyForHook, void * > hook;

	// if 'this' is a face of another cell and that other cell belongs to some mesh msh,
	// cell_behind_within[msh] keeps that cell
	// see methods Mesh::cell_behind and Mesh::cell_in_front_of
	std::map < Mesh::Core *, Cell > cell_behind_within;
	// we use the Cell wrapper as pointer

	// we use a wrapper as pointer
	Cell reverse_attr;

	inline CellCore ( const tag::OfDimension &, const size_t d, const tag::OneDummyWrapper & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Core::DelegateDispose ( tag::one_dummy_wrapper ),
	#else  // no MANIFEM_COLLECT_CM
	:
	#endif  // MANIFEM_COLLECT_CM	
		double_heap ( Cell::double_heap_size_pos [d] ),
		size_t_heap ( Cell::size_t_heap_size_pos [d] ),
		short_int_heap ( Cell::short_int_heap_size_pos [d] )
	{ }

	inline CellCore ( const tag::OfDimension &, const size_t d, const tag::ZeroWrappers & )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Core::DelegateDispose ( tag::zero_wrappers ),
	#else  // no MANIFEM_COLLECT_CM
	:
	#endif  // MANIFEM_COLLECT_CM	
		double_heap ( Cell::double_heap_size_neg [d] ),
		size_t_heap ( Cell::size_t_heap_size_neg [d] ),
		short_int_heap ( Cell::short_int_heap_size_neg [d] )
	{ }

	virtual ~CellCore ( ) { };

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	virtual bool is_positive ( ) const = 0;
	virtual Cell get_positive ( ) = 0;
	virtual size_t get_dim ( ) const = 0;

	virtual Cell tip ();
	virtual Cell base ();
	// execution forbidden
	// overridden by Cell::Positive::Segment and Cell::Negative::Segment

	virtual Mesh boundary ( ) = 0;

	virtual bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const = 0;
	virtual bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const = 0;

	// Method 'glue_on_bdry_of' is intensively used when building a mesh,
	// e.g. within factory functions in Cell class.
	// It glues 'this' cell to the boundary of 'cll'.
	inline void glue_on_bdry_of ( Cell::Core * cll )
	{	cll->glue_on_my_bdry ( this );   }

	// Method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll'. Used mainly in remeshing.
	inline void cut_from_bdry_of ( Cell::Core * cll )
	{	cll->cut_from_my_bdry ( this );   }
	
	// the two methods below are only relevant for vertices
	// so we forbid execution for now and then override them in Cell::***tive::Vertex
  virtual void add_to_seg ( Cell::PositiveSegment * seg );
	virtual void remove_from_seg ( Cell::PositiveSegment * seg );

	// the four methods below are not relevant for vertices
  virtual void add_to_mesh ( Mesh::Core * msh ) = 0;
	virtual void remove_from_mesh ( Mesh::Core * msh ) = 0;
  virtual void add_to_bdry ( Mesh::Core * msh ) = 0;
	virtual void remove_from_bdry ( Mesh::Core * msh ) = 0;

	virtual void glue_on_my_bdry ( Cell::Core * ) = 0;
	virtual void cut_from_my_bdry ( Cell::Core * ) = 0;
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see above)

	#ifndef NDEBUG
	std::string name;
	virtual std::string get_name ( ) = 0;
	virtual void print_everything ( ) = 0;
	#endif

}; // end of class Cell::Core

//-----------------------------------------------------------------------------//


class Cell::Positive : public Cell::Core

// abstract class, introduces attribute  meshes
// and defines virtual methods  is_positive, get_positive, belongs_to

// specialized in Cell::Positive{Vertex,Segment,HighDim}

// I would very much prefer the name
// class 'Cell::Positive::Vertex' instead of 'Cell::PositiveVertex'
// but it was not possible ...

{	public :

	typedef Cell::PositiveVertex Vertex;
	typedef Cell::PositiveNotVertex NotVertex;
	typedef Cell::PositiveSegment Segment;
	typedef Cell::PositiveHighDim HighDim;
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// the 'meshes' attribute keeps information about all meshes
	// "above" 'this' cell, that is, all meshes containing 'this' as a cell
	// it is indexed over the dimension of the mesh minus the dimension of 'this' 
	// for each mesh, it keeps a 'Cell::field_to_meshes' value, containing two counters
	// and an iterator into the 'cells' field of that particular mesh
	// the iterator is only meaningful for fuzzy or STSI meshes
	// of course this implies quite some amount of redundant information
	// but this redundancy makes the classes fast, especially for remeshing
	
	// since indices of vectors begin at zero
	// the keys of 'meshes[i]' will be meshes of dimension 'i + this->dim'
	// however, for meshes of the same dimension as 'this' cell,
	// we keep a different record in Cell::Positive::Vertex::segments
	// and Cell::Positive::NotVertex::meshes_same_dim
	// so meshes[0] will always be an empty map

	std::vector < std::map < Mesh::Core*, Cell::field_to_meshes > > meshes;

	inline Positive ( const tag::OfDimension &, const size_t d,
	                  const tag::SizeMeshes &, const size_t sz, const tag::OneDummyWrapper & )
	:	Cell::Core ( tag::of_dim, d, tag::one_dummy_wrapper ),
		meshes ( sz )
	{	}

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell get_positive ( );  // virtual from Cell::Core

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	#ifndef NDEBUG
	std::string get_name();
	#endif
	
}; // end of class Cell::Positive


class Cell::Negative : public Cell::Core

// abstract class, defines virtual methods
// is_positive, get_positive, build_reverse, belongs_to,
// glue_on_my_bdry, cut_from_my_bdry

// specialized in Cell::Negative::{Vertex,Segment,HighDim}

{	public :

	typedef Cell::NegativeVertex Vertex;
	typedef Cell::NegativeSegment Segment;
	typedef Cell::NegativeHighDim HighDim;
	// I would very much prefer the name 'Cell::Negative::Vertex', 'Cell::Negative::Segment'
	// instead of 'Cell::NegativeVertex', 'Cell::NegativeSegment'
	// but it was not possible ...
	
	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	inline Negative ( const tag::OfDimension &, const size_t d,
	                  const tag::ReverseOf &, Cell::Positive * rev, const tag::OneDummyWrapper & )
	: Cell::Core ( tag::of_dim, d, tag::one_dummy_wrapper ),
		reverse_attr ( tag::whose_core_is, rev, tag::previously_existing )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell get_positive ( );  // virtual from Cell::Core
	// virtual size_t get_dim ( ) const = 0;  // declared in Cell::Core

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

	#ifndef NDEBUG
	std::string get_name();
	#endif
	
}; // end of class Cell::Negative

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::Vertex' instead of 'Cell::PositiveVertex'
// but it was not possible ...

class Cell::PositiveVertex : public Cell::Positive

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell::Negative * reverse_p  inherited from Cell::Positive

	// meshes  inherited from Cell::Positive

	// in 'segments' we keep record of meshes of dimension zero "above" 'this' vertex,
	// that is, of segments having 'this' extremity
	// the 'short int' is a sign, 1 or -1

	std::map < Cell::Positive::Segment *, short int > segments;
	
	inline PositiveVertex ( const tag::OneDummyWrapper & );

	PositiveVertex ( const Cell::Positive::Vertex & ) = delete;
	PositiveVertex ( const Cell::Positive::Vertex && ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex & ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex && ) = delete;

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	// two versions of  belongs_to  defined by Cell::Positive

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// the two methods below are defined by Cell::Core but overriden here
	void add_to_seg ( Cell::Positive::Segment * seg ) override;
	void remove_from_seg ( Cell::Positive::Segment * seg ) override;

	// the six methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh );
	void glue_on_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core * );

	#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
	#endif
	
}; // end of  class Cell::PositiveVertex, aka Cell::Positive::Vertex


// I would very much prefer the name
// 'class Cell::Negative::Vertex' instead of 'Cell::NegativeVertex'
// but it was not possible ...

class Cell::NegativeVertex : public Cell::Negative

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Negative

	inline NegativeVertex
	( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p, const tag::OneDummyWrapper & );

	NegativeVertex ( const Cell::Negative::Vertex & ) = delete;
	NegativeVertex ( const Cell::Negative::Vertex && ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex & ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex && ) = delete;

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core but overriden here
	void add_to_seg ( Cell::Positive::Segment * seg ) override;
	void remove_from_seg ( Cell::Positive::Segment * seg ) override;

	// the four methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh );

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::NegativeVertex, aka Cell::Negative::Vertex

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::NotVertex' instead of 'Cell::PositiveNotVertex'
// but it was not possible ...

class Cell::PositiveNotVertex : public Cell::Positive

// abstract class, useful only for introducing the attribute  meshes_same_dim
// and methods glue_common, cut_common
// specialized in Cell::Positive::{Segment,HighDim}

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// meshes  inherited from Cell::Positive

	// in 'meshes_same_dim' we keep record of meshes "above" 'this' cell,
	//   of the same dimension
	// in Cell::field_to_meshes_same_dim, the 'short int sign' is a sign, 1 or -1
	// the iterator 'where' is only meaningful for fuzzy or STSI meshes

	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > meshes_same_dim;
	
	inline PositiveNotVertex ( const tag::OfDimension &, const size_t d,
                             const tag::SizeMeshes &, const size_t sz,
                             const tag::OneDummyWrapper &              )
	: Cell::Positive ( tag::of_dim, d, tag::size_meshes, sz, tag::one_dummy_wrapper )
	{	}
	
	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	inline void glue_common ( Cell::Core * face );
	inline void cut_common ( Cell::Core * face );
	// do not use directly; called from glue_on_my_bdry and cut_from_my_bdry

}; // end of  class Cell::PositiveNotVertex, aka Cell::Positive::NotVertex

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::Segment' instead of 'Cell::PositiveSegment'
// but it was not possible ...

class Cell::PositiveSegment : public Cell::Positive::NotVertex

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell::Negative * reverse_p  inherited from Cell::Positive

	// meshes  inherited from Cell::Positive
	// meshes_same_dim  inherited from Cell::Positive::NotVertex

	// here we use Cell wrappers as pointers
	// because a segment should keep its extremities alive
	Cell base_attr, tip_attr;

	inline PositiveSegment
		( Cell::Negative::Vertex * a, Cell::Positive::Vertex * b, const tag::OneDummyWrapper & );

	PositiveSegment ( const Cell::Positive::Segment & ) = delete;
	PositiveSegment ( const Cell::Positive::Segment && ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment & ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	Cell tip () override; // virtual, overrides definition by Cell::Core
	Cell base () override; // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

	// glue_common  and  cut_common  defined by Cell::Positive::NotVertex

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::PositiveSegment, aka Cell::Positive::Segment


// I would very much prefer the name
// 'class Cell::Negative::Segment' instead of 'Cell::NegativeSegment'
// but it was not possible ...

class Cell::NegativeSegment : public Cell::Negative

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Negative

	inline NegativeSegment
	( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p, const tag::OneDummyWrapper & );

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	NegativeSegment ( const Cell::Negative::Segment & ) = delete;
	NegativeSegment ( const Cell::Negative::Segment && ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment & ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	Cell tip () override;  // virtual, overrides definition by Cell::Core
	Cell base () override;  // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::NegativeSegment, aka Cell::Negative::Segment

//-----------------------------------------------------------------------------//


// I would very much prefer the name
// 'class Cell::Positive::HighDim' instead of 'Cell::PositiveHighDim'
// but it was not possible ...

class Cell::PositiveHighDim : public Cell::Positive::NotVertex

// a cell of dimension >= 2

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell::Negative * reverse_p  inherited from Cell::Positive

	// meshes  inherited from Cell::Positive
	// meshes_same_dim  inherited from Cell::Positive::NotVertex

	// we use Mesh wrapper as pointer
	// because a Cell should keep its boundary alive
	Mesh boundary_attr;  // positive mesh

	inline PositiveHighDim ( const tag::OfDimension &, const size_t d,
	                         const tag::WhoseBoundaryIs &, Mesh msh, const tag::OneDummyWrapper & );
	inline PositiveHighDim ( const tag::WhoseBoundaryIs &, Mesh msh, const tag::OneDummyWrapper & );
	inline PositiveHighDim ( const tag::Triangle &, Cell::Core * AB,
	                         Cell::Core * BC, Cell::Core * CA, const tag::OneDummyWrapper &        );
	inline PositiveHighDim ( const tag::Quadrangle &, Cell::Core * AB, Cell::Core * BC,
	                         Cell::Core * CD, Cell::Core * DA, const tag::OneDummyWrapper & );

	PositiveHighDim ( const Cell::Positive::HighDim & ) = delete;
	PositiveHighDim ( const Cell::Positive::HighDim && ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim & ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim && ) = delete;

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual

	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Positive

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

	// glue_common  and  cut_common  defined by Cell::Positive::NotVertex
	
#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

}; // end of  class Cell::PositiveHighDim, aka Cell::Positive::HighDim


// I would very much prefer the name
// 'class Cell::Negative::HighDim' instead of 'Cell::NegativeHighDim'
// but it was not possible ...

class Cell::NegativeHighDim : public Cell::Negative

// a negative cell of dimension >= 2

{	public :

	// size_t nb_of_wrappers  inherited from tag::Util::Core ifdef COLLECT_CM

	// double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// Cell reverse_attr  inherited from Cell::Negative

	inline NegativeHighDim ( const tag::OfDimension &, const size_t d,
	                         const tag::ReverseOf &, Cell::Positive::HighDim * direct_cell_p, const tag::OneDummyWrapper & );
	
	inline NegativeHighDim ( const tag::ReverseOf &, Cell::Positive::HighDim * direct_cell_p, const tag::OneDummyWrapper & );

	NegativeHighDim ( const Cell::Negative::HighDim & ) = delete;
	NegativeHighDim ( const Cell::Negative::HighDim && ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim & ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim && ) = delete;

	// bool dispose ( )  defined by tag::Util::Core::DelegateDispose ifdef COLLECT_CM
	
	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core

	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

}; // end of  class Cell::NegativeHighDim, aka Cell::Negative::HighDim


//-----------------------------------------------------------------------------//
//------------------       core Meshes     ------------------------------------//
//-----------------------------------------------------------------------------//


// I would very much prefer the name 'class Mesh::Core' instead of 'tag::Util::MeshCore'
// but it was not possible ...

#ifdef MANIFEM_COLLECT_CM	

class tag::Util::MeshCore : public tag::Util::Core

#else  // no MANIFEM_COLLECT_CM

class tag::Util::MeshCore

#endif  // MANIFEM_COLLECT_CM	

// represents a positive mesh
// negative meshes have no core (wrappers for negative meshes are built on-the-fly)

// abstract class, specialized in Mesh::ZeroDim, Mesh::Connected::{OneDim,HighDim},
// Mesh::MultiplyConnected::{OneDim,HighDim}, Mesh::Fuzzy, Mesh::STSI

{	public :

	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// we use an ordinary pointer, not a wrapper
	// because we do not want a cell to be kept alive by its boundary
	Cell::Positive * cell_enclosed { nullptr };
	// when 'this' is the boundary of a cell

	inline MeshCore ( const tag::OfDimension &, const size_t d , const tag::MinusOne &,
                    const tag::OneDummyWrapper &                                      )
	#ifdef MANIFEM_COLLECT_CM	
	:	tag::Util::Core ( tag::one_dummy_wrapper )
	#endif  // MANIFEM_COLLECT_CM	
	{	}
	
	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	virtual size_t get_dim_plus_one ( ) = 0;

	virtual size_t number_of ( const tag::Vertices & ) = 0;
	virtual size_t number_of ( const tag::Segments & ) = 0;
	virtual size_t number_of ( const tag::CellsOfDim &, const size_t d ) = 0;
	virtual size_t number_of ( const tag::CellsOfMaxDim & ) = 0;

	// the four methods below are only relevant for connected one-dimensional meshes
	// so we forbid execution for now and then override them in Mesh::Connected::OneDim
	virtual Cell first_vertex ( );
	virtual Cell last_vertex ( );
	virtual Cell first_segment ( );
	virtual Cell last_segment ( );

	// the four methods below are only relevant for STSI meshes
	// so we forbid execution for now and then override them in Mesh::STSI
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                               ) const;
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::MayNotExist &                                                       ) const;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                               ) const;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::MayNotExist &                                                       ) const;

	virtual void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;

	virtual Mesh::Core * build_deep_copy ( ) = 0;
	
	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core *, const size_t );
	// returns garbage; overriden by Mesh::Fuzzy and later by Mesh::STSI

	virtual void remove_from_cells
	( Cell::Core *, const size_t, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// does nothing; overriden by Mesh::Fuzzy and later by Mesh::STSI
		
	// iterators defined in iterator.cpp
	// we are still in class Mesh::Core

	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   ) = 0;

	// we are still in class Mesh::Core

	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   ) = 0;

	// we are still in class Mesh::Core

	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive &               ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;

	// we are still in class Mesh::Core

	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & ) = 0;

	// several versions of 'build' below are defined in global.cpp

	void build ( const tag::Segment &,  // builds a chain of n segment cells
	             const Cell & A, const Cell & B, const tag::DividedIn &, const size_t n );

	void build ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	
	void build ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	
	void build ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	             const Cell & NE, const Cell & NW, const size_t m, const size_t n,
	             bool cut_rectangles_in_half                                 );
	
#ifndef NDEBUG
	std::string name;
  virtual std::string get_name() = 0;
 	virtual void print_everything () = 0;
#endif

}; // end of  class Mesh::Core


class Mesh::ZeroDim : public Mesh::Core

// represents two points, one negative one positive

// zero-dimensional meshes only exist as boundary of a segment,
// which can be retrieved through the cell_enclosed attribute
// they are built on-the-fly

{	public :

	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed inherited from Mesh::Core
	
	inline ZeroDim ( const tag::OneDummyWrapper & )
	:	Mesh::Core ( tag::of_dimension, 1, tag::minus_one, tag::one_dummy_wrapper )
	{ }

	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );
	// virtual from Mesh::Core, here execution forbidden
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private:

	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// the sixteen methods below are virtual from Mesh::Core, here execution forbidden
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	
	Mesh::Core * build_deep_copy ( );  // virtual from Mesh::Core, here execution forbidden

	// we are still in class Mesh::ZeroDim
	
	// add_to_cells ( Cell::Core *, const size_t )  defined by Cell::Core, returns garbage

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base (negative) then tip (positive)

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first tip (positive) then base (negative)

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	// iterate over the two vertices, first base then tip (both positive)

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first tip then base (both positive)

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::ZeroDim

	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::ZeroDim
	// twelve iterators below assert dim == 0 then call iterator over vertices

	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive &               );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::ZeroDim
	// twelve iterators below simply produce iterator over vertices
	
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

};  // end of class Mesh::ZeroDim


class Mesh::Connected::OneDim : public Mesh::Core

// represents a connected positive mesh of dimension 1
// a chain of segments, either open or closed

{	public :

	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed inherited from Mesh::Core

	size_t nb_of_segs;
	// useful for quickly answering to 'number_of'

	// here we use Cell wrappers as pointers
	Cell first_ver, last_ver;
	// if last_ver == first_ver->reverse, the mesh is a closed loop
	
	// inline OneDim ( )
	//	:	Mesh::Core ( tag::of_dimension, 2, tag::minus_one ) { }

	OneDim ( const tag::Segment &, Cell::Negative::Vertex * A, Cell::Positive::Vertex * B,
	            const tag::DividedIn &, const size_t n, const tag::OneDummyWrapper & );
	// defined in global.cpp
	
	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	Cell first_vertex ( );  // virtual from Mesh::Core, here overridden
	Cell last_vertex ( );  // virtual from Mesh::Core, here overridden
	Cell first_segment ( );  // virtual from Mesh::Core, here overridden
	Cell last_segment ( );  // virtual from Mesh::Core, here overridden

	// private:

	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, const size_t )  defined by Cell::Core, returns garbage

	inline void order ( );
	// run over all segments, order them linearly, check that the mesh is connected
	// if it is a loop, set first_ver = Cell::ghost
	// if it is an open chain, set first_ver and last_ver
	// see paragraph 9.14 in the manual
	
#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

}; // end of  class Mesh::Connected::OneDim


class Mesh::Connected::HighDim : public Mesh::Core

// represents a connected positive mesh of dimension >= 2

{	public :

	// attribute cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// Cell start;
	
	inline HighDim ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &, const tag::OneDummyWrapper & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

  HighDim ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA, const tag::OneDummyWrapper & );
	// defined in global.cpp

	inline HighDim ( const tag::Quadrangle &,
		const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
									 ,const tag::OneDummyWrapper &,
		const tag::WithTriangles & wt = tag::not_with_triangles                       );

	inline HighDim ( const tag::Quadrangle &,
		const Cell & SW, const Cell & SE, const Cell & NE, const Cell & NW,
		const size_t m, const size_t n, const tag::OneDummyWrapper &,
									 const tag::WithTriangles & wt = tag::not_with_triangles );

	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, const size_t )  defined by Cell::Core, returns garbage

	inline CellIterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::Connected::HighDim


class Mesh::MultiplyConnected::OneDim : public Mesh::Core

// represents a positive mesh of dimension 1 with several connected components

{	public :

	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting/ending points

	inline OneDim ( const tag::OneDummyWrapper & )
		:	Mesh::Core ( tag::of_dimension, 2, tag::minus_one, tag::one_dummy_wrapper )
	{ }

  OneDim ( const tag::Segment &, Cell::Negative::Vertex * A,
		Cell::Positive::Vertex * B, const tag::DividedIn &, const size_t n, const tag::OneDummyWrapper & );
	// defined in global.cpp
	
	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	// private:
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, const size_t )  defined by Cell::Core, returns garbage

	inline void order ( );
	// run over all segments, order them linearly, check that the mesh is connected
	// if it is a loop, set first_ver = Cell::ghost
	// if it is an open chain, set first_ver and last_ver
	// see paragraph 9.14 in the manual
	
#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

}; // end of  class Mesh::MultiplyConnected::OneDim


class Mesh::MultiplyConnected::HighDim : public Mesh::Core

// represents a positive mesh of dimension >= 2 with several connected components

{	public :

	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting cells
	
	inline HighDim ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &, const tag::OneDummyWrapper & )
		:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	static bool is_positive ( );
	static Mesh reverse ( Mesh::Core * core );

	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, const size_t )  defined by Cell::Core, returns garbage

	inline CellIterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::MultiplyConnected::HighDim


class Mesh::Fuzzy : public Mesh::Core

// represents a positive mesh, disordered (includes 1D meshes)

// roughly speaking, a mesh is a collection of cells of the same dimension
// however, for efficiency purposes, we keep lists of cells of lower
// dimension as well; that's the purpose of the 'cells' vector :
// to keep lists of cells indexed by their dimension

{	public :
	
	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed inherited from Mesh::Core

	// the 'cells' attribute holds lists of cells of 'this' mesh, indexed by their dimension
	// for maximum dimension, the cells are oriented
	// for lower dimension, the cells are always positive

	// here we use Cell wrappers as pointers
	// because a mesh should keep its cells alive
	std::vector < std::list < Cell > > cells;

	inline Fuzzy ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &, const tag::OneDummyWrapper & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper ), cells ( dim_p1 )
	{	assert ( dim_p1 > 1 );  }
		
	inline Fuzzy ( const tag::Segment &,  // builds a chain of n segment cells
	               const Cell & A, const Cell & B, const tag::DividedIn &, const size_t n, const tag::OneDummyWrapper & )
	:	Fuzzy ( tag::of_dimension, 2, tag::minus_one, tag::one_dummy_wrapper )
	{	build ( tag::segment, A, B, tag::divided_in, n );  }
	
	inline Fuzzy ( const tag::Triangle &,  // builds a triangle from its three sides
	               const Mesh & AB, const Mesh & BC, const Mesh & CA, const tag::OneDummyWrapper & )
	:	Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper )
	{	build ( tag::triangle, AB, BC, CA );  }

	inline Fuzzy ( const tag::Quadrangle &,  // builds a rectangle from its four sides
	               const Mesh & south, const Mesh & east,
                 const Mesh & north, const Mesh & west, const tag::OneDummyWrapper &,
	               const tag::WithTriangles & wt = tag::not_with_triangles )
	// 'wt' defaults to 'tag::not_with_triangles'
	:	Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper )
	{	build ( tag::quadrangle, south, east, north, west, wt == tag::with_triangles );  }
	
	inline Fuzzy ( const tag::Quadrangle &,  // builds a rectangle from its four sides
	               const Cell & SW, const Cell & SE, const Cell & NE, const Cell & NW,
	               const size_t m, const size_t n, const tag::OneDummyWrapper &,
	               const tag::WithTriangles & wt = tag::not_with_triangles              );
	
	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, const size_t d );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfMaxDim & );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp
	
	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// the sixteen methods below are virtual from Mesh::Core
	// called from Cell::****tive::***::add_to_mesh and Cell::****tive::***::remove_from_mesh
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );

	Mesh::Core * build_deep_copy ( );  // virtual from Mesh::Core

	// add a cell to 'this->cells[d]' list, return iterator into that list
	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core * cll, const size_t d ) override;
	// virtual from Cell::Core, here overriden, later overriden again by Mesh::STSI
	
	// remove a cell from 'this->cells[d]' list using the provided iterator
	virtual void remove_from_cells
	( Cell::Core *, const size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// virtual from Cell::Core, here overriden, later overriden again by Mesh::STSI
	
	// we are still in class Mesh::Fuzzy
	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverVertices &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::ZeroDim

	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &                        );
	CellIterator::Core * iterator
	( const tag::OverSegments &, const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &                   );

	// we are still in class Mesh::Fuzzy

	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &        );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &        );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &   );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive &            );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &            );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &       );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive &               );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfDim &, const size_t,
	  const tag::ReverseEachCell &, const tag::DoNotBother &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );

	// we are still in class Mesh::Fuzzy
	
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
	  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator::Core * iterator
	( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
	  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::Fuzzy


class Mesh::STSI : public Mesh::Fuzzy

// represents a positive mesh, maybe self-touching, maybe self-intersecting
// (includes 1D meshes)

// mainly used for iterators over connected meshes
// also used for progressive mesh generation

{	public :
	
	// unsigned int nb_of_wrappers  inherited from tag::Util::Core

	// cell_enclosed  inherited from Mesh::Core

	// 'cells' inherited from Mesh::Fuzzy

	// in 'singular' we keep pair of adjacent cells
	// the common face of such a pair is a singular face
	// that is, a face where the mesh touches itself
	std::vector < std::pair < Cell, Cell > > singular;
	
	inline STSI ( const tag::OfDimension &, const size_t dim_p1, const tag::MinusOne &, const tag::OneDummyWrapper & )
	:	Mesh::Fuzzy ( tag::of_dimension, dim_p1, tag::minus_one, tag::one_dummy_wrapper )
	{	}

	// bool dispose ( )  defined by tag::Util::Core ifdef COLLECT_CM
	
	// size_t get_dim_plus_one ( )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Vertices & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Segments & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::CellsOfDim &, const size_t d )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::CellsOfMaxDim & )  defined in Mesh::Fuzzy

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp
	
	static Mesh reverse ( Mesh::Core * core );

	// private :

	// the four methods below are defined in Cell::Core, here overriden
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &,
    const Cell::Core * neighbour, const tag::MayNotExist & ) const override;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &,
    const Cell::Core * neighbour, const tag::MayNotExist & ) const override;

	// add a cell to 'this->cells[d]' list, return iterator into that list
	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core * cll, const size_t d ) override;
	// virtual from Cell::Core, overriden by Mesh::Fuzzy, here overriden again

	// remove a cell from 'this->cells[d]' list using the provided iterator
	virtual void remove_from_cells
	( Cell::Core *, const size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// virtual from Cell::Core, overriden by Mesh::Fuzzy, here overriden again

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	inline CellIterator iterator ( const tag::OverCellsOfDim &, const size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::STSI

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


inline void Mesh::copy_all_cells_to ( Mesh & msh ) const
{	CellIterator it = this->iterator ( tag::over_cells_of_max_dim );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell cll = *it;  cll.add_to ( msh );  }                         }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                    const tag::IsNegative &, const tag::CellsSurelyExist &                     )
// builds a negative mesh from a positive one, assuming that reverse cells exist
// used in Cell::boundary and in Mesh::Mesh below

: Mesh ( tag::whose_core_is, c, tag::previously_existing,
         tag::is_negative, tag::do_not_bother             )

#ifndef NDEBUG
{	assert ( c );
	// check that all cells have reverse
	CellIterator it = this->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell::Core * cll_p = (*it).core;
		Cell::Positive * cll_pos = dynamic_cast < Cell::Positive* > ( cll_p );
		if ( cll_pos )
			assert ( cll_pos->reverse_p );
		else
		{ Cell::Negative * cll_neg = dynamic_cast < Cell::Negative* > ( cll_p );
			assert ( cll_neg );
			assert ( cll_neg->reverse_attr.exists() );
			assert ( cll_neg->reverse_attr.reverse ( tag::surely_exists ) == *it );    }  }  }
#else
{	}
#endif


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::PreviouslyExisting &,
                    const tag::IsNegative &, const tag::BuildCellsIfNec & b                    )
// builds a negative mesh from a positive one, creating reverse cells if necessary
// used in Mesh::Positive::reverse
	
: Mesh ( tag::whose_core_is, c, tag::previously_existing,
         tag::is_negative, tag::do_not_bother             )

{	CellIterator it = this->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset() ; it.in_range(); it++ )
		it->reverse ( tag::build_if_necessary );            }

inline Mesh::Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t d,
                    const tag::IsPositive & ispos                                 )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	Mesh ( tag::whose_core_is,
	 	     new Mesh::Fuzzy ( tag::of_dimension, d+1, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                      )
{	}


inline Mesh::Mesh ( const tag::Fuzzy &, const tag::OfDimension &, const size_t d,
                    const tag::MinusOne &, const tag::IsPositive & ispos          )
// by default, ispos = tag::is_positive, so may be called with only four arguments
:	Mesh ( tag::whose_core_is,
	 	     new Mesh::Fuzzy ( tag::of_dimension, d, tag::minus_one, tag::one_dummy_wrapper ),
	       tag::freshly_created, tag::is_positive                      )
{	}
	

inline Mesh::Mesh ( const tag::DeepCopyOf &, const Mesh & msh )
:	Mesh ( tag::whose_core_is, msh.core->build_deep_copy(),
	       tag::freshly_created, tag::is_positive           )
// this may be dangerous, the core is built with tag::one_dummy_wrapper
// we must ensure that it does not meet any other wrapper before 'this' one
// in other words, we must ensure that the core is freshly created
// and eager to meet its parent (wrapper)
{	}

inline Mesh::Mesh ( const tag::DeepCopyOf &, const Mesh & msh, const tag::Fuzzy & )
:	Mesh ( tag::fuzzy, tag::of_dim, msh.core->get_dim_plus_one(),
         tag::minus_one, tag::is_positive                       )
{	msh.copy_all_cells_to ( *this );  }


inline Mesh::Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
                    const tag::DividedIn &, const size_t n                      )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::segment, A, B, tag::divided_in, n, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                      )
{	}


inline Mesh::Mesh
(	const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::triangle, AB, BC, CA, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive         )
{	}


inline Mesh::Mesh ( const tag::Quadrangle &, const Mesh & south,
                    const Mesh & east, const Mesh & north, const Mesh & west,
                    const tag::WithTriangles & wt                             )
// 'wt' defaults to 'tag::not_with_triangles',
// so constructor can be called with only five arguments
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::quadrangle, south, east, north, west, tag::one_dummy_wrapper, wt ),
         tag::freshly_created, tag::is_positive                             )
{	}

inline Mesh::Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
                    const Cell & NE, const Cell & NW, const size_t m, const size_t n,
                    const tag::WithTriangles & wt                                     )
// 'wt' defaults to 'tag::not_with_triangles',
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::quadrangle, SW, SE, NE, NW, m ,n, tag::one_dummy_wrapper, wt ),
         tag::freshly_created, tag::is_positive                         )
{	}


inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2 )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, m1.core->get_dim_plus_one(), tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                        )
{	std::vector < Mesh > l;  l.reserve ( 2 );
	l.push_back ( m1 );  l.push_back ( m2 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3 )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, m1.core->get_dim_plus_one(), tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                        )
{	std::vector < Mesh > l;  l.reserve ( 3 );
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
  this->join_list ( l );                                         }

inline Mesh::Mesh
( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3, const Mesh & m4 )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, m1.core->get_dim_plus_one(), tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                        )
{	std::vector < Mesh > l;  l.reserve ( 4 );
	l.push_back ( m1 );  l.push_back ( m2 );
	l.push_back ( m3 );  l.push_back ( m4 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3,
                                       const Mesh & m4, const Mesh & m5 )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, m1.core->get_dim_plus_one(), tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                        )
{	std::vector < Mesh > l;  l.reserve ( 5 );
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
	l.push_back ( m4 );  l.push_back ( m5 );
  this->join_list ( l );                                         }

template < typename container >
inline Mesh::Mesh ( const tag::Join &, const container & l )
: Mesh ( tag::whose_core_is,
	       new Mesh::Fuzzy ( tag::of_dim, l.front().core->get_dim_plus_one(), tag::minus_one, tag::one_dummy_wrapper ),
         tag::freshly_created, tag::is_positive                                               )
{	this->join_list ( l );  }


//-----------------------------------------------------------------------------//


inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh & msh )
:	Cell ( tag::whose_core_is,
	       new Cell::Positive::HighDim ( tag::whose_boundary_is, msh, tag::one_dummy_wrapper ),
				 tag::freshly_created                                         )
#ifndef NDEBUG
{	assert ( msh.is_positive() );
	assert ( msh.dim() >= 1 );
	Mesh::STSI * msh_stsi = dynamic_cast < Mesh::STSI* > ( msh.core );
	assert ( msh_stsi == nullptr );                                     }
#else
{	}
#endif

inline Cell::Cell ( const tag::Vertex &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex, tag::freshly_created )
{	}


inline Cell::Cell ( const tag::Segment &, const Cell & A, const Cell & B )
: Cell ( tag::whose_core_is, new Cell::Positive::Segment
	        ( tag::Util::assert_cast < Cell::Core*, Cell::NegativeVertex* > ( A.core ),
	          tag::Util::assert_cast < Cell::Core*, Cell::PositiveVertex* > ( B.core ),
						tag::one_dummy_wrapper                                                    ),
	       tag::freshly_created                                                           )
{	}


inline Cell::Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
				 ( tag::triangle, AB.core, BC.core, CA.core, tag::one_dummy_wrapper ), tag::freshly_created )
{	}

inline Cell::Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                             const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive::HighDim
				 ( tag::quadrangle, AB.core, BC.core, CD.core, DA.core, tag::one_dummy_wrapper ), tag::freshly_created )
{	}

//-----------------------------------------------------------------------------//


inline size_t Mesh::number_of ( const tag::Vertices & ) const
{	return this->number_of ( tag::cells_of_dim, 0 );  }

inline size_t Mesh::number_of ( const tag::Segments & ) const
{	return this->number_of ( tag::cells_of_dim, 1 );  }

inline size_t Mesh::number_of ( const tag::CellsOfDim &, const size_t d ) const
{	return this->core->number_of ( tag::cells_of_dim, d );  }

inline size_t Mesh::number_of ( const tag::CellsOfMaxDim & ) const
{	return this->core->number_of ( tag::cells_of_max_dim );  }


inline Cell Mesh::first_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_vertex().core,
		              tag::previously_existing                             );
	// else
	return Cell ( tag::whose_core_is,
	              this->core->last_vertex().core->reverse ( tag::surely_exists ),
                tag::previously_existing                                        );  }

inline Cell Mesh::last_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_vertex().core,
		              tag::previously_existing                             );
	// else
	return Cell ( tag::whose_core_is,
	              this->core->first_vertex().core->reverse ( tag::surely_exists ),
                tag::previously_existing                                         );  }

inline Cell Mesh::first_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_segment().core,
		              tag::previously_existing                              );
	// else
	return Cell ( tag::whose_core_is,
	              this->core->last_vertex().core->reverse ( tag::surely_exists ),
                tag::previously_existing                                        );  }

inline Cell Mesh::last_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_segment().core,
		              tag::previously_existing                             );
	// else
	return Cell ( tag::whose_core_is,
	              this->core->first_vertex().core->reverse ( tag::surely_exists ),
                tag::previously_existing                                         );  }


inline bool Cell::belongs_to ( const Mesh & msh, const tag::Oriented & ) const
// what about negative meshes ? haa ?!!	
{	return this->core->belongs_to ( msh.core, tag::oriented );  }

inline bool Cell::belongs_to ( const Mesh & msh, const tag::NotOriented & ) const
// what about negative meshes ? haa ?!!	
{	return this->core->belongs_to ( msh.core, tag::not_oriented );  }

//-----------------------------------------------------------------------------//


inline Mesh::Fuzzy::Fuzzy ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
														const Cell & NE, const Cell & NW, const size_t m, const size_t n, tag::one_dummy_wrapper, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

	:	Mesh::Fuzzy ( tag::of_dimension, 3, tag::minus_one, tag::one_dummy_wrapper )
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	Mesh south ( tag::segment, SW.reverse(), SE, tag::divided_in, m );
	Mesh east  ( tag::segment, SE.reverse(), NE, tag::divided_in, n );
	Mesh north ( tag::segment, NE.reverse(), NW, tag::divided_in, m );
	Mesh west  ( tag::segment, NW.reverse(), SW, tag::divided_in, n );

	this->build ( tag::rectangle, south, east, north, west, cut_rectangles_in_half );  }

//-----------------------------------------------------------------------------//


inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell face_rev = face.reverse ( tag::may_not_exist );
	if ( not face_rev.exists() ) return nullptr;
	else return this->cell_behind ( face_rev, tag::may_not_exist );  }
	

inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell face_rev = face.reverse ( tag::surely_exists );
	assert ( face_rev.exists() );
	return this->cell_behind ( face_rev, tag::surely_exists );      }
	

inline Cell Mesh::cell_behind
( const Cell face, const tag::MayNotExist & ) const

// return the cell to which 'face' belongs, non-existent if we are facing the boundary

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell>::const_iterator
			it = face.core->cell_behind_within.find ( this->core );
		if ( it == face.core->cell_behind_within.end() ) return Cell ( tag::non_existent );
			// nothing behind us, we are touching the boundary
		assert ( it->second.exists() );
		return it->second;                                                                  }
		// face.core->cell_behind_within[this->core]
	else
	{	Cell face_rev = face.reverse ( tag::surely_exists );
		// we are in a negative mesh, all faces must have reverse
		std::map<Mesh::Core*,Cell>::const_iterator
			it = face_rev.core->cell_behind_within.find ( this->core );
		if ( it == face_rev.core->cell_behind_within.end() )
			return nullptr;  // we are facing the boundary
		Cell cll_rev = it->second;
		assert ( cll_rev.exists() );
		return cll_rev.reverse ( tag::surely_exists );                  }  }


inline Cell Mesh::cell_behind ( const Cell face, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell to which 'face' belongs

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell>::const_iterator
			it = face.core->cell_behind_within.find ( this->core );
		assert ( it != face.core->cell_behind_within.end() );
		assert ( it->second.exists() );
		return it->second;                                        }
		// face.core->cell_behind_within[this->core]
	else
	{	Cell face_rev = face.reverse ( tag::surely_exists );
		std::map<Mesh::Core*,Cell>::const_iterator
			it = face_rev.core->cell_behind_within.find ( this->core );
		assert ( it != face_rev.core->cell_behind_within.end() );
		Cell cll_rev = it->second;
		assert ( cll_rev.exists );
		return cll_rev.reverse ( tag::surely_exists );                 }  }


// the four methods below are only relevant for STSI meshes

inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::SeenFrom &,
  const Cell::Core * neighbour, const tag::SurelyExists & se ) const
// return the cell towards which 'face' points, as seen from neighbour
// 'se' defaults to tag::surely_exists, so method can be called with only three arguments
{	return this->core->cell_in_front_of ( face_p, tag::seen_from, neighbour, se );  }

inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
  const tag::MayNotExist &                                                        ) const
// return the cell towards which 'face' points, as seen from neighbour
{	return this->core->cell_in_front_of
		( face_p, tag::seen_from, neighbour, tag::may_not_exist );  }

inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SeenFrom &,
  const Cell::Core * neighbour, const tag::SurelyExists & se ) const
// return the cell to which 'face' belongs, as seen from neighbour
// 'se' defaults to tag::surely_exists, so method can be called with only three arguments
{	return this->core->cell_behind ( face_p, tag::seen_from, neighbour, se );  }

inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core * neighbour,
  const tag::MayNotExist &                                                        ) const
// return the cell to which 'face' belongs, as seen from neighbour
{	return this->core->cell_behind ( face_p, tag::seen_from, neighbour, tag::may_not_exist );  }


#ifndef NDEBUG
inline void Mesh::print_everything ( )
{	if ( not is_positive() ) std::cout << "(negative Mesh) ";
	core->print_everything ();                                }
#endif

//-----------------------------------------------------------------------------//


inline bool Cell::is_positive ( ) const
{	return this->core->is_positive ( );  }

inline Cell Cell::get_positive ( )
{	return Cell ( tag::whose_core_is, this->core->get_positive(), tag::previously_existing );  }


inline size_t Cell::dim ( ) const
{	return this->core->get_dim ( );  }

inline size_t Mesh::dim ( ) const
{	return tag::Util::assert_diff ( this->core->get_dim_plus_one(), 1 );      }
// tag::Util::assert_diff  provides a safe way to substract two size_t numbers


inline Cell Cell::reverse ( const tag::BuildIfNotExists & build ) const
// 'build' defaults to tag::build_if_not_exists, so method may be called with no arguments
{	if ( not this->core->reverse_attr.exists() )
		this->core->reverse_attr =
			Cell ( tag::whose_core_is, this->core->build_reverse(), tag::freshly_created )
	return this->core->reverse_attr;                                                  }

inline Cell Cell::reverse ( const tag::MayNotExist & ) const
{	return this->core->reverse_attr;  }

inline Cell Cell::reverse ( const tag::SurelyExists & ) const
{	assert ( this->core->reverse_attr.exists() );
	return this->core->reverse_attr;              }

inline Mesh Mesh::reverse ( ) const
{	if ( this->is_positive() )
		return Mesh ( tag::whose_core_is, core, tag::previously_existing,
	                tag::is_negative, tag::build_cells_if_necessary     );
	// else
	return Mesh ( tag::whose_core_is, core, tag::previously_existing, tag::is_positive );  }
	

inline bool Cell::has_reverse ( ) const
{	return this->core->reverse ( tag::may_not_exist );  }


inline Mesh Cell::boundary ( ) const

{	assert ( this->core );
	return this->core->boundary();    }

inline Cell Cell::tip () const
{	// assert ( this->core->tip() );
	return Cell ( tag::whose_core_is, this->core->tip().core, tag::previously_existing );  }

inline Cell Cell::base () const
{	// assert ( this->core->base() );
	return Cell ( tag::whose_core_is, this->core->base().core, tag::previously_existing );  }


#ifndef NDEBUG
inline void Cell::print_everything ( )
{	core->print_everything ( );  }
#endif


inline void Cell::glue_on_bdry_of ( Cell & cll )

// glue 'this' face on the boundary of cell 'cll'
// any of them may be negative

{	this->core->glue_on_bdry_of ( cll.core );  }


inline void Cell::cut_from_bdry_of ( Cell & cll )

// cut 'this' face from the boundary of cell 'cll'
// any of them may be negative

{	this->core->cut_from_bdry_of ( cll.core );  }


inline void Cell::add_to ( Mesh & msh )

// add 'this' cell to the mesh 'msh' by calling the virtual method add_to_mesh
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( this->dim() == msh.dim() );
	assert ( this->dim() > 0 );
	if ( msh.is_positive() )  this->core->add_to_mesh ( msh.core );
	else
		this->core->reverse ( tag::surely_exists ) ->add_to_mesh ( msh.core );  }
// for negative Meshes, the core points towards the reverse Mesh::Core


inline void Cell::remove_from ( Mesh & msh )

// remove 'this' cell from the mesh 'msh' by calling the virtual method remove_from_mesh
// if 'msh' is the boundary of a cell, 'cut_from_bdry_of' should be used instead

{	assert ( this->dim() == msh.dim() );
	assert ( this->dim() > 0 );
	if ( msh.is_positive() )  this->core->remove_from_mesh ( msh.core );
	else
		this->core->reverse ( tag::surely_exists ) ->remove_from_mesh ( msh.core );  }
// for negative Meshes, the core points towards the reverse Mesh::Core


//-----------------------------------------------------------------------------//


inline Cell::PositiveVertex::PositiveVertex ( const tag::OneDummyWrapper & )
	: Cell::Positive ( tag::of_dim, 0, tag::size_meshes, Mesh::maximum_dimension_plus_one, tag::one_dummy_wrapper )
{	}


inline Cell::NegativeVertex::NegativeVertex
( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p, const tag::OneDummyWrapper & )
	: Cell::Negative ( tag::of_dim, 0, tag::reverse_of, direct_ver_p, tag::one_dummy_wrapper )
{	}


inline Cell::PositiveSegment::PositiveSegment
( Cell::Negative::Vertex * Aa, Cell::Positive::Vertex * Bb, const tag::OneDummyWrapper & )

: Cell::Positive::NotVertex ( tag::of_dim, 1, tag::size_meshes,
		  tag::Util::assert_diff ( Mesh::maximum_dimension_plus_one, 1 ) ),
	// tag::Util::assert_diff provides a safe way to substract two 'size_t' numbers
	base_attr ( tag::whose_core_is, Aa, tag::previously_existing ),
	tip_attr  ( tag::whose_core_is, Bb, tag::previously_existing )

{	// below is a much simplified version of Cell::Negative::Vertex::add_to_seg
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the base has already been correctly initialized
	Cell::Positive::Vertex * pos_Aa = tag::Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( Aa->reverse ( tag::surely_exists ) );
	assert ( pos_Aa->segments.find(this) == pos_Aa->segments.end() );
	// pos_Aa->meshes[0][msh] = Cell::field_to_meshes { 0, 1 };
	// the third component 'where' is irrelevant here
	pos_Aa->segments.emplace ( std::piecewise_construct,
	      std::forward_as_tuple(this), std::forward_as_tuple(-1) );
	// below is a much simplified version of Cell::Positive::Vertex::add_to_seg
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the tip has already been correctly initialized
	assert ( Bb->segments.find(this) == Bb->segments.end() );
	// Bb->meshes[0][msh] = Cell::field_to_meshes { 1, 0 };
	// the third component 'where' is irrelevant here
	Bb->segments.emplace ( std::piecewise_construct,
	   std::forward_as_tuple(this), std::forward_as_tuple(1) );                         }


inline Cell::NegativeSegment::NegativeSegment
( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p, const tag::OneDummyWrapper & )

: Cell::Negative ( tag::of_dim, 1, tag::reverse_of, direct_seg_p, tag::one_dummy_wrapper )

// we must make sure that both extremities of 'direct_seg_p' have a reverse
// well, the base surely has one since it's a NegativeVertex

{	assert ( direct_seg_p->base_attr.exists() );
	assert ( direct_seg_p->tip_attr.exists() );
	assert ( direct_seg_p->base_attr.core->reverse ( tag::surely_exists ) );
	Cell::Core * seg_rev_p = direct_seg_p->tip_attr.core->reverse ( tag::build_if_not_exists );
	assert ( seg_rev_p == direct_seg_p->tip_attr.core->reverse ( tag::surely_exists ) );        }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::OfDimension &, const size_t d,
  const tag::WhoseBoundaryIs &, Mesh msh, const tag::OneDummyWrapper &    )

:	Cell::Positive::NotVertex ( tag::of_dim, d, tag::size_meshes,
	         tag::Util::assert_diff ( Mesh::maximum_dimension_plus_one, d ), tag::one_dummy_wrapper ),
	// tag::Util::assert_diff provides a safe way to substract two 'size_t' numbers
	boundary_attr ( msh )
	
{	assert ( msh.core );
	assert ( msh.core->get_dim_plus_one() == d );
	msh.core->cell_enclosed = this;           }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::WhoseBoundaryIs &, Mesh msh, const tag::OneDummyWrapper & )
:	Cell::Positive::HighDim ( tag::of_dim, msh.core->get_dim_plus_one(),
	                          tag::whose_bdry_is, msh, tag::one_dummy_wrapper )
{	}


inline Cell::NegativeHighDim::NegativeHighDim
( const tag::ReverseOf &, Cell::Positive::HighDim * direct_cell_p, const tag::OneDummyWrapper & )
:	Cell::Negative::HighDim
	( tag::of_dim, direct_cell_p->get_dim(), tag::reverse_of, direct_cell_p, tag::one_dummy_wrapper )
{	}


inline Cell::NegativeHighDim::NegativeHighDim
(	const tag::OfDimension &, const size_t d, const tag::ReverseOf &,
	Cell::PositiveHighDim * direct_cell_p, const tag::OneDummyWrapper & )
	
: Cell::Negative ( tag::of_dim, d, tag::reverse_of, direct_cell_p, tag::one_dummy_wrapper )

// we must make sure that all faces of 'direct_cell_p' have a reverse  #ifndef NDEBUG !!

{	assert ( direct_cell_p );
	assert ( direct_cell_p->get_dim() == d );
	assert ( direct_cell_p->boundary.exists() );
	assert ( direct_cell_p->boundary.is_positive() );
	CellIterator it = direct_cell_p->boundary.iterator
		( tag::over_cells_of_max_dim, tag::as_they_are );
	// CellIterator it ( tag::whose_core_is, direct_cell_p->boundary_p->iterator
	// 	( tag::over_cells_of_max_dim, tag::as_they_are, tag::this_mesh_is_positive )  );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * cll_p = (*it).core;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p == cll_p->reverse ( tag::surely_exists ) );        }         }


inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Triangle &, Cell::Core * AB, Cell::Core * BC, Cell::Core * CA, const tag::OneDummyWrapper & )
	
: Cell::Positive::HighDim ( tag::whose_boundary_is,
														new Mesh::Fuzzy ( tag::of_dim, 2, tag::minus_one, tag::one_dummy_wrapper ), tag::freshly_created )
// in the future, boundary will be Mesh::Connected::OneDim !!
	
{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CA->base()->reverse_p );
	assert ( CA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to_bdry ( this->boundary_p );
	BC->add_to_bdry ( this->boundary_p );
	CA->add_to_bdry ( this->boundary_p );            }
	

inline Cell::PositiveHighDim::PositiveHighDim
( const tag::Quadrangle &, Cell::Core * AB, Cell::Core * BC, Cell::Core * CD, Cell::Core * DA, const tag::OneDummyWrapper & )
	
: Cell::Positive::HighDim ( tag::whose_boundary_is,
														new Mesh::Fuzzy ( tag::of_dim, 2, tag::minus_one, tag::one_dummy>wrapper ), tag::freshly_created )
// in the future, boundary will be Mesh::Connected::OneDim !!

{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CD->get_dim() == 1 );
	assert ( DA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CD->base()->reverse_p );
	assert ( CD->tip() == DA->base()->reverse_p );
	assert ( DA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to_bdry ( this->boundary_p );
	BC->add_to_bdry ( this->boundary_p );
	CD->add_to_bdry ( this->boundary_p );
	DA->add_to_bdry ( this->boundary_p );            }

//---------------------------------------------------------------------------------


inline void Cell::Positive::NotVertex::glue_common ( Cell::Core * face )

{	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > & tm0 = this->meshes_same_dim;
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim>::iterator it;
	for ( it = tm0.begin(); it != tm0.end(); ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
			face->cell_behind_within[msh] = this;
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			rev_face->cell_behind_within[msh] = other_cell;   }     }                            }
	

inline void Cell::Positive::NotVertex::cut_common ( Cell::Core * face )
	
{	if ( this->meshes.size() == 0 ) return;
	std::map < Mesh::Core *, Cell::field_to_meshes_same_dim > & tm0 = this->meshes_same_dim;
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim>::iterator it;
	for ( it = tm0.begin(); it != tm0.end(); ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
		{	assert ( face->cell_behind_within.find(msh) !=
		           face->cell_behind_within.end()        );
			assert ( face->cell_behind_within[msh] == this );
			face->cell_behind_within.erase(msh);                  }
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			assert ( rev_face->cell_behind_within.find(msh) !=
		           rev_face->cell_behind_within.end()        );
			assert ( rev_face->cell_behind_within[msh] == other_cell );
			rev_face->cell_behind_within.erase(msh);                     }  }                  }
	


}  // namespace maniFEM

#endif
// ifndef MANIFEM_MESH_H
