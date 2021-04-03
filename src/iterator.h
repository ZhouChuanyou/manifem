
// iterator.h 2021.04.02

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

#ifndef MANIFEM_ITERATOR_H
#define MANIFEM_ITERATOR_H

#include "mesh.h"

namespace maniFEM {

namespace tag {
	struct OverTwoVerticesOfSeg { };
	static const OverTwoVerticesOfSeg over_two_vertices_of_seg;
	struct OverVerticesOf { };  static const OverVerticesOf over_vertices_of;
	struct OverSegmentsOf { };  static const OverSegmentsOf over_segments_of;
	struct OverCellsOf { };  static const OverCellsOf over_cells_of;
	struct FuzzyPosMesh { };  static const FuzzyPosMesh fuzzy_pos_mesh;
	struct ReverseEachCell { };  static const ReverseEachCell reverse_each_cell;
	struct AssumeCellsExist { };  static const AssumeCellsExist assume_cells_exist;
}


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
	
	inline CellIterator ( const tag::OverTwoVerticesOfSeg &,
	                      Cell::PositiveSegment * seg         );
	// we iterate over the two vertices, first base (negative) then tip (positive)
		
	inline CellIterator ( const tag::OverTwoVerticesOfSeg &,
	                      Cell::PositiveSegment * seg, const tag::ReverseOrder & );
	// we iterate over the two vertices, first tip (positive) then base (negative)

	inline CellIterator ( const tag::OverTwoVerticesOfSeg &,
	                      Cell::PositiveSegment * seg, const tag::ForcePositive & );
	// we iterate over the two vertices, first base then tip (both positive)

	inline CellIterator ( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
	                      const tag::ReverseOrder &, const tag::ForcePositive &           );
	// we iterate over the two vertices, first tip then base (both positive)
	
	inline CellIterator ( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
	                      const tag::ReverseEachCell &, const tag::AssumeCellsExist &     );
	// we iterate over the two vertices, first base (positive) then tip (negative)

	inline CellIterator ( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
                        const tag::ReverseOrder &,
	                      const tag::ReverseEachCell &, const tag::AssumeCellsExist &     );
	// we iterate over the two vertices, first tip (negative) then base (positive)

	inline CellIterator ( const tag::OverCellsOf &, Mesh::Fuzzy * msh,
	                      const tag::FuzzyPosMesh &,
	                      const tag::CellsOfDim &, size_t, const tag::AsTheyAre & );

	inline CellIterator ( const tag::OverCellsOf &, Mesh::Fuzzy * msh,
	                      const tag::FuzzyPosMesh &,
	                      const tag::CellsOfDim &, size_t d, const tag::ForcePositive &  );

	inline CellIterator ( const tag::OverCellsOf &, Mesh::Fuzzy * msh,
	                      const tag::FuzzyPosMesh &,
	                      const tag::CellsOfMaxDim &, const tag::AsTheyAre & );

	inline CellIterator ( const tag::OverCellsOf &, Mesh::Fuzzy * msh,
	                      const tag::FuzzyPosMesh &,
	                      const tag::CellsOfMaxDim &, const tag::ForcePositive & );

	inline CellIterator ( const tag::OverCellsOf &, Mesh::Fuzzy * msh,
	                      const tag::FuzzyPosMesh &,
	                      const tag::CellsOfMaxDim &, const tag::ReverseEachCell & );

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


class CellIterator::Core

// iterates over all cells of a given mesh (cells of given dimension)
// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative
// for lower dimension, returns positive cells

// or, iterates over all cells above a given cell

{	public :

	inline Core ( ) { };
	
	virtual ~Core ( ) { };

	virtual void reset ( ) = 0;
	virtual void reset ( Cell::Core * cll ) = 0;
	virtual Cell::Core * deref ( ) = 0;
	virtual void advance ( ) = 0;
	virtual bool in_range ( ) = 0;
	
};  // end of class CellIterator::Core


inline CellIterator::CellIterator ( const CellIterator & it)
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy construct an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                }
// define virtual self-replicator for core iterators ?
	
inline CellIterator::CellIterator ( CellIterator && it)
{	this->core = std::move ( it.core );  } // it.~CellIterator();  }
// should we destroy 'it' ?  if we do,
// the core gets destroyed because that's how we define CellIterator::~CellIterator below
// then how does the new iterator work ? without a core ?
// if we don't, what happens to the memory slot occupied by 'it', more precisely by it->core ?
// is it automatically released ?  hope so ...
	
inline CellIterator & CellIterator::operator= ( const CellIterator & it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to copy an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// define virtual self-replicator for core iterators ?

inline CellIterator & CellIterator::operator= ( const CellIterator && it )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "ManiFEM doesn't know (yet) how to move an iterator. Sorry." << std::endl;
	exit ( 1 );                                                                                   }
// ownership transfer ? I don't know how to implement it ...
// come on bro, it's easy, just copy the pointer (or rather, move it)
	
inline void CellIterator::reset ( )  {  this->core->reset();  }

inline void CellIterator::reset ( Cell & cll )
{  this->core->reset( cll.core );  }
	
inline Cell CellIterator::operator* ( )
{	return Cell ( tag::whose_core_is, this->core->deref() );  }

inline CellIterator & CellIterator::operator++ ( )  {  return this->advance();  }

inline CellIterator & CellIterator::operator++ ( int )  {  return this->advance();  }

inline CellIterator & CellIterator::advance ( )
{ this->core->advance();  return *this;  }

inline bool CellIterator::in_range ( )  {  return this->core->in_range();  }
		

class CellIterator::Over::TwoVerticesOfSeg : public CellIterator::Core

{	public :

	Cell::PositiveSegment * seg_p;
	size_t passage;

	inline TwoVerticesOfSeg ( Cell::PositiveSegment * seg )
	:	CellIterator::Core (), seg_p { seg } { };
	
	void reset ( );  // virtual from CellIterator::Core
	void reset ( Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden
	// Cell::Core * deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	struct NormalOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	struct ReverseOrder  {  class AsTheyAre;  class ForcePositive;
		struct ReverseEachCell  {  class  AssumeCellsExist;  };    };
	
};  // end of class CellIterator::Over::TwoVerticesOfSeg


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

	class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline AsTheyAre ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

	class ForcePositive;
	struct ReverseEachCell  {  class AssumeCellsExist;  };

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline ForcePositive ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive


class CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::NormalOrder::ReverseEachCell::AssumeCellsExist


class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ReverseEachCell::AssumeCellsExist
: public CellIterator::Over::TwoVerticesOfSeg

{	public :

	inline AssumeCellsExist ( Cell::PositiveSegment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg
    //         ::ReverseOrder::ReverseEachCell::AssumeCellsExist


class CellIterator::Over::CellsOfFuzzyMesh : public CellIterator::Core

{	public :

	std::list<Cell::Core*> * list_p;
  std::list<Cell::Core*>::iterator iter;

	inline CellsOfFuzzyMesh ( Mesh::Fuzzy * msh, size_t dim )
	:	CellIterator::Core (), list_p { & ( msh->cells[dim] ) }
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	void reset ( );  // virtual from CellIterator::Core
	void reset ( Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden
	// Cell::Core * deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class AsTheyAre;  class ForcePositive;  class ReverseEachCell;
	
};  // end of class CellIterator::Over::CellsOfFuzzyMesh


class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	std::list<Cell::Core*> * list_p;
  std::list<Cell::Core*>::iterator iter;

	inline AsTheyAre ( Mesh::Fuzzy * msh, size_t dim )
	:	CellIterator::Over::CellsOfFuzzyMesh ( msh, dim )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre


class CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	std::list<Cell::Core*> * list_p;
  std::list<Cell::Core*>::iterator iter;

	inline ReverseEachCell ( Mesh::Fuzzy * msh, size_t dim )
	:	CellIterator::Over::CellsOfFuzzyMesh ( msh, dim )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// we trust each cell has already a reverse
	// void advance ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell


class CellIterator::Over::CellsOfFuzzyMesh::ForcePositive
: public CellIterator::Over::CellsOfFuzzyMesh

{	public :

	std::list<Cell::Core*> * list_p;
  std::list<Cell::Core*>::iterator iter;

	inline ForcePositive ( Mesh::Fuzzy * msh, size_t dim )
	:	CellIterator::Over::CellsOfFuzzyMesh ( msh, dim )
	{	}  // no need to initialize 'iter', 'reset' will do that
	
	// void reset ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::CellsOfFuzzyMesh, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::CellsOfFuzzyMesh
	// bool in_range ( )  defined by CellIterator::Over::CellsOfFuzzyMesh

};  // end of class CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre

//--------------------------------------------------------------------------------

	
class CellIterator::Adaptor::ForcePositive : public CellIterator::Core

// modifies the behaviour of another CellIterator
// forcing the resulting cells to be positive

{	public :

	CellIterator base;

	inline ForcePositive ( CellIterator && b )
	:	base { b } { };

	void reset ( );  // virtual from CellIterator::Core
	void reset ( Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden
	// Cell::Core * deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

};  // end of class CellIterator::Adaptor::ForcePositive

//-------------------------------------------------------------------------------------------


inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg )
// we iterate over the two vertices, first base (negative) then tip (positive)
:	core { new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg ) }
{ }
		
inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg, const tag::ReverseOrder & )
// we iterate over the two vertices, first tip (positive) then base (negative)
:	core { new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg ) }
{ }

inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg, const tag::ForcePositive & )
// we iterate over the two vertices, first base then tip (both positive)
:	core { new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg ) }
{ }

inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
  const tag::ReverseOrder &, const tag::ForcePositive &              )
// we iterate over the two vertices, first tip then base (both positive)
:	core { new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg ) }
{ }
		
inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
  const tag::ReverseEachCell &, const tag::AssumeCellsExist &      )
// we iterate over the two vertices, first base (positive) then tip (negative)
:	core { new CellIterator::Over::TwoVerticesOfSeg
		::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg ) }
{ }
		
inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfSeg &, Cell::PositiveSegment * seg,
  const tag::ReverseOrder &, const tag::ReverseEachCell &, const tag::AssumeCellsExist & )
// we iterate over the two vertices, first tip (negative) then base (positive)
:	core { new CellIterator::Over::TwoVerticesOfSeg
		::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg ) }
{ }
		
inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfDim &, size_t d, const tag::AsTheyAre &                 )
: core { new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( d )  }
{ }

inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfDim &, size_t d, const tag::ReverseEachCell &           )
: core { new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell ( d )  }
{ }

inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfDim &, size_t d, const tag::ForcePositive &           )
: core { new CellIterator::Over::CellsOfFuzzyMesh::ForcePositive ( d )  }
{ }

inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfMaxDim &, const tag::AsTheyAre &                      )
: core { new CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre ( d )  }
{ }

inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfMaxDim &, const tag::ForcePositive &                  )
: core { new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell ( d )  }
{ }

inline CellIterator::CellIterator
( const tag::OverCellsOf &, Mesh::Fuzzy * msh, const tag::FuzzyPosMesh &,
  const tag::CellsOfMaxDim &, const tag::ReverseEachCell &                )
: core { new CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell ( d )  }
{ }

//-------------------------------------------------------------------------------------------


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::CellsSurelyExist &                                      )
// builds a negative mesh from a positive one, assuming that reverse cells exist
// used in Cell::boundary and in Mesh::Mesh below

: Mesh ( tag::whose_core_is, c, tag::is_negative, tag::cells_may_not_exist )

#ifndef NDEBUG
{	assert ( c );
	// check that all cells have reverse
	CellIterator it = this->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell::Core * cll_p = (*it).core;
		assert ( cll_p );  assert ( cll_p->reverse_p );
		assert ( cll == cll_p->reverse_p->reverse_p );   }           }
#else
{	}
#endif


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::BuildCellsIfNec & b                                     )
// builds a negative mesh from a positive one, creating reverse cells if necessary
// used in Mesh::Positive::reverse
	
:	Mesh ( tag::whose_core_is, c, tag::is_negative, tag::cells_may_not_exist )

{	CellIterator it = this->iterator  // as they are : oriented
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell::Core * cll_p = (*it).core;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p );  assert ( cll_rev_p == cll_p->reverse_p );             }     }


inline Cell::NegativeHighDim::NegativeHighDim
( const tag::OfDimension &, size_t d, const tag::ReverseOf &, Cell::PositiveHighDim * direct_cell_p )
	
: Cell::Negative ( tag::of_dim, d, tag::reverse_of, direct_cell_p )

// we must make sure that all faces of 'direct_cell_p' have a reverse

{	assert ( direct_cell_p );
	assert ( direct_cell_p->get_dim() == d );
	assert ( direct_cell_p->boundary_p );
	CellIterator it = direct_cell_p->boundary_p->iterator
		( tag::over_cells_of_max_dim, tag::as_they_are );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * cll_p = (*it).core;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p == cll_p->reverse_p );  assert ( cll_rev_p );             }    }


//-------------------------------------------------------------------------------------------


inline CellIterator Mesh::iterator ( const tag::OverVertices & ) const
{	if ( this->is_positive() )
		return this->core->iterator ( tag::over_vertices, tag::this_mesh_is_positive );
	// else
	return this->core->iterator ( tag::over_vertices, tag::this_mesh_is_negative );    }

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
		( tag::over_vertices, tag::force_positive, tag::this_mesh_is_positive );    }

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseOrder & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_vertices, tag::reverse_order, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
		( tag::over_vertices, tag::reverse_order, tag::this_mesh_is_negative );     }

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ReverseOrder & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_vertice, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
			( tag::over_vertices, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_negative );  }

inline CellIterator Mesh::iterator
( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ForcePositive & ) const
{	return this->iterator
		( tag::over_vertices, tag::force_positive, tag::reverse_order );  }

inline CellIterator Mesh::iterator ( const tag::OverSegments & ) const
{	if ( this->is_positive() )
		return this->core->iterator ( tag::over_segments, tag::this_mesh_is_positive );
	// else
	return this->core->iterator ( tag::over_segments, tag::this_mesh_is_negative );    }

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_segments, tag::force_positive, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
		( tag::over_segments, tag::force_positive, tag::this_mesh_is_positive );    }

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseOrder & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_segments, tag::reverse_order, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
		( tag::over_segments, tag::reverse_order, tag::this_mesh_is_negative );     }

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ReverseOrder & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_segments, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
			( tag::over_segments, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_negative );  }

inline CellIterator Mesh::iterator
( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ForcePositive & ) const
{	return this->iterator
		( tag::over_segments, tag::force_positive, tag::reverse_order );  }

inline CellIterator Mesh::iterator	( const tag::OverCellsOfDim &, size_t d ) const
{	if ( this->is_positive() )
		return this->core->iterator ( tag::over_cells_of_dim, d, tag::this_mesh_is_positive );
	// else
	return this->core->iterator ( tag::over_cells_of_dim, d, tag::this_mesh_is_negative );   }

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder & ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_cells_of_dim, d, tag::reverse_order, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
		( tag::over_cells_of_dim, d, tag::reverse_order, tag::this_mesh_is_negative );    }

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive & ) const
{	if ( d < this->dim() )
		return this->core->iterator ( tag::over_cells_of_max_dim, tag::force_positive );
	return this->core->iterator ( tag::over_cells_of_dim, d, tag::lower_dimension );   }

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
  const tag::ReverseOrder &                                           ) const
{	if ( this->is_positive() )
		return this->core->iterator
			( tag::over_cells_of_dim, d, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_positive );
	// else
	return this->core->iterator
			( tag::over_cells_of_dim, d, tag::reverse_order,
			  tag::force_positive, tag::this_mesh_is_negative );  }

inline CellIterator Mesh::iterator
( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
  const tag::ForcePositive &                                         ) const
{	return this->iterator
		( tag::over_cells_of_dim, d, tag::force_positive, tag::reverse_order );  }

//-----------------------------------------------------------------------------//


inline CellIterator Mesh::Core::iterator
( const tag::OverVertices, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return this->iterator ( tag::over_vertices, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_positive );  }

inline CellIterator Mesh::Core::iterator
( const tag::OverVertices, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsNegative & )
{	return this->iterator ( tag::over_vertices, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_negative );  }

inline CellIterator Mesh::Core::iterator
( const tag::OverSegments, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	return this->iterator ( tag::over_segments, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_positive );  }

inline CellIterator Mesh::Core::iterator
( const tag::OverSegments, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsNegative & )
{	return this->iterator ( tag::over_segments, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_negative );  }

inline CellIterator Mesh::Core::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsPositive &               )
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_positive );  }

inline CellIterator Mesh::Core::iterator
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
	const tag::ReverseOrder &, const tag::ThisMeshIsNegative &               )
{	return this->iterator ( tag::over_cells_of_dim, d, tag::reverse_order,
										      tag::force_positive, tag::this_mesh_is_negative );  }

}  // namespace maniFEM

#endif
// ifndef MANIFEM_ITERATOR_H
