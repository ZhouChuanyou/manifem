
// iterator.h 2021.03.23

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
	struct OverTwoVerticesOfPosSeg;
	static const OverTwoVerticesOfPosSeg over_two_vertices_of_pos_seg;
	struct OverTwoVerticesOfNegSeg;
	static const OverTwoVerticesOfNegSeg over_two_vertices_of_neg_seg;
}


class CellIterator

// a thin wrapper around a CellIterator::Core, with most methods delegated to 'core'

// iterates over all cells of a given mesh (cells of given dimension)

// for maximum dimension (equal to the dimension of the mesh) returns oriented cells
// which may be positive or negative
// for lower dimension, returns positive cells

// or, iterates over all cells above a given cell

{	public :

	class Core;
	
	std::unique_ptr < CellIterator::Core > core;

	inline CellIterator ( const CellIterator & it );
	inline CellIterator ( CellIterator && it );

	inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &,
	                      const tag::ForcePositive &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &,
	                      const tag::ReverseOrder &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &,
    const tag::ReverseOrder &, const tag::ForcePositive &,, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &,
	                      const tag::ForcePositive &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &,
	                      const tag::ReverseOrder &, Cell::Positive::Segment * );
	inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &,
    const tag::ReverseOrder &, const tag::ForcePositive &,, Cell::Positive::Segment * );

	inline ~CellIterator ( ) { };
	
	inline void reset ( );	
	inline void reset ( Cell & );	
	inline Cell operator* ( );
	inline CellIterator & operator++ ( );
	inline CellIterator & operator++ ( int );
	inline CellIterator & advance ( );
	inline bool in_range ( );

	struct Over
	{	class TwoVerticesOfSeg;  };
	struct TwoVerticesOfPosSeg  {  class NormalOrder;  class ReverseOrder;  };
	struct TwoVerticesOfNegSeg  {  class NormalOrder;  class ReverseOrder;  };
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
{	this->core = it.core;  } // it.~CellIterator();  }
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
{  this->core = this->core->reset( cll.core );  }
	
inline Cell CellIterator::operator* ( )
{	return Cell ( tag::whose_core_is, this->core->deref() );  }

inline CellIterator & CellIterator::operator++ ( )  {  return this->advance();  }

inline CellIterator & CellIterator::operator++ ( int )  {  return this->advance();  }

inline CellIterator & CellIterator::advance ( )
{ this->core->advance();  return *this;  }

inline bool CellIterator::in_range ( )  {  return this->core->in_range();  }
		

class CellIterator::Over::TwoVerticesOfSeg : public CellIterator::Core

{	public :

	Cell::Positive::Segment * seg_p;
	size_t try;

	inline TwoVerticesOfSeg ( Cell::Positive::Segment * seg )
	:	seg_p { seg } { };
	
	void reset ( );  // virtual from CellIterator::Core
	void reset ( Cell::Core * cll );
	// virtual from CellIterator::Core, here execution forbidden
	// Cell::Core * deref ( )  remains pure virtual from CellIterator::Core
	void advance ( );  // virtual from CellIterator::Core
	bool in_range ( );  // virtual from CellIterator::Core

	class NormalOrder;  class ReverseOrder;
	
};  // end of class CellIterator::Over::TwoVerticesOfSeg


class CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder
: public CellIterator::TwoVerticesOfSeg

{	public :

	inline NormalOrder ( Cell::Positive::Segment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

	class ForcePositive;

};  // end of class CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder


class CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder::ForcePositive
: public CellIterator::TwoVerticesOfSeg

{	public :

	inline ForcePositive ( Cell::Positive::Segment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::NormalOrder


class CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder
: public CellIterator::TwoVerticesOfSeg

{	public :

	inline ReverseOrder ( Cell::Positive::Segment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

	class ForcePositive;

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder


class CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder::ForcePositive
: public CellIterator::TwoVerticesOfSeg

{	public :

	inline ReverseOrder ( Cell::Positive::Segment * seg )
	:	CellIterator::Over::TwoVerticesOfSeg ( seg ) { };
	
	// void reset ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// void reset ( Cell::Core * cll )
	// defined by CellIterator::Over::TwoVerticesOfSeg, execution forbidden
	Cell::Core * deref ( );  // virtual from CellIterator::Core
	// void advance ( )  defined by CellIterator::Over::TwoVerticesOfSeg
	// bool in_range ( )  defined by CellIterator::Over::TwoVerticesOfSeg

};  // end of class CellIterator::Over::TwoVerticesOfSeg::ReverseOrder


inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfPosSeg &, Cell::Positive::Segment * seg )
// we iterate over the two vertices, first base (negative) then tip (positive)
:	core { new CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder ( seg ) }
{ }
	
inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &,
                      const tag::ReverseOrder &, Cell::Positive::Segment * )
// we iterate over the two vertices, first tip (positive) then base (negative)
:	core { new CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder ( seg ) }
{ }
	
inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &,
                      const tag::ForcePositive &, Cell::Positive::Segment * )
// we iterate over the two vertices, first base then tip (both positive)
:	core { new CellIterator::Over::TwoVerticesOfPosSeg::ForcePositive ( seg ) }
{ }

inline CellIterator ( const tag::OverTwoVerticesOfPosSeg &, const tag::ReverseOrder &,
                      const tag::ForcePositive &, Cell::Positive::Segment *            )
// we iterate over the two vertices, first tip then base (both positive)
:	core { new CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder::ForcePositive ( seg ) }
{ }

inline CellIterator::CellIterator
( const tag::OverTwoVerticesOfNegSeg &, Cell::Positive::Segment * seg )
// we iterate over the two vertices, first base (negative) then tip (positive)
:	core { new CellIterator::Over::TwoVerticesOfNegSeg::NormalOrder ( seg ) }
{ }
	
inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &,
                      const tag::ReverseOrder &, Cell::Positive::Segment * )
// we iterate over the two vertices, first tip (positive) then base (negative)
:	core { new CellIterator::Over::TwoVerticesOfNegSeg::ReverseOrder ( seg ) }
{ }
	
inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &,
                      const tag::ForcePositive &, Cell::Positive::Segment * )
// we iterate over the two vertices, first base then tip (both positive)
:	core { new CellIterator::Over::TwoVerticesOfNegSeg::ForcePositive ( seg ) }
{ }
	
inline CellIterator ( const tag::OverTwoVerticesOfNegSeg &, const tag::ReverseOrder &,
                      const tag::ForcePositive &, Cell::Positive::Segment *            )
// we iterate over the two vertices, first tip then base (both positive)
:	core { new CellIterator::Over::TwoVerticesOfNegSeg::ReverseOrder::ForcePositive ( seg ) }
{ }


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

}  // end of class CellIterator::Adaptor::ForcePositive


}  // namespace maniFEM

#endif
// ifndef MANIFEM_ITERATOR_H
