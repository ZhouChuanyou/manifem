
// iterator.cpp 2021.04.05

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

#include "iterator.h"

using namespace maniFEM;

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );     }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );     }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                 }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                  }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                      }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                              }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Zero-dimensional meshes have no segments." << std::endl;
	exit ( 1 );                                                             }


CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }
	
CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );      }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &          )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );      }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive &              )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive &              )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );   }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive &         )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );   }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ThisMeshIsPositive &                 )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                 }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                  }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( d == 0 );
	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (negative) then tip (positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::AsTheyAre ( seg );       }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );     }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::AsTheyAre &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (positive) then base (negative)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::AsTheyAre ( seg );     }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base then tip (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::NormalOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ForcePositive &,
  const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip then base (both positive)
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg::ReverseOrder::ForcePositive ( seg );    }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                 }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::RequireOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first base (positive) then tip (negative)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::NormalOrder::ReverseEachCell::AssumeCellsExist ( seg );                  }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }

CellIterator::Core * Mesh::ZeroDim::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfMaxDim &, const tag::ReverseEachCell &,
  const tag::DoNotBother &, const tag::ReverseOrderIfAny &, const tag::ThisMeshIsPositive & )
// a positive zero-dimensional mesh is the boundary of a positive segment	
// iterate over the two vertices, first tip (negative) then base (positive)
// do not bother whether reverse cells exist or not
{	assert ( this->cell_enclosed );
	Cell::Positive::Segment * seg =
		assert_cast < Cell::Positive *, Cell::Positive::Segment * > ( this->cell_enclosed );
	return new CellIterator::Over::TwoVerticesOfSeg
	             ::ReverseOrder::ReverseEachCell::AssumeCellsExist ( seg );                }


CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ThisMeshIsPositive & )
{	return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
                        tag::cells_of_dim, 0, tag::as_they_are         );   }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
// fuzzy meshes have dimension at least one
// thus, vertices are positive anyway
{	return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
                        tag::cells_of_dim, 0, tag::as_they_are         );   }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverVertices &, const tag::ReverseOrder &,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > 1  );
	// fuzzy meshes have dimension at least one
	return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
												tag::cells_of_dim, 1, tag::as_they_are          );  }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > 1 );
	// fuzzy meshes have dimension at least one
	if ( this->get_dim_plus_one() == 2 )  // one-dimensional fuzzy mesh
		return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
													tag::cells_of_dim, 1, tag::force_positive       );
	else  // high-dimensional fuzzy mesh
		return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
													tag::cells_of_dim, 1, tag::as_they_are          );  }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverSegments &, const tag::ReverseOrder &,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive & )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ThisMeshIsPositive & )
{	assert ( this->get_dim_plus_one() > d );
	return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
												tag::cells_of_dim, 1, tag::as_they_are          );  }
	
CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseOrder &,
  const tag::ThisMeshIsPositive &                                   )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ForcePositive &,
  const tag::ThisMeshIsPositive &                                          )
{	assert ( this->get_dim_plus_one() > d );
	if ( this->get_dim_plus_one() == d+1 )  // cells of maximum dimension
		return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
													tag::cells_of_dim, d, tag::force_positive       );
	else  // cells of lower dimension
		return CellIterator ( tag::over_cells_of, this, tag::fuzzy_pos_mesh,
													tag::cells_of_dim, d, tag::as_they_are          );  }

CellIterator::Core * Mesh::Fuzzy::iterator  // virtual from Mesh::Core
( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseOrder &,
  const tag::ForcePositive &, const tag::ThisMeshIsPositive &             )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "It makes no sense to reverse order for a fuzzy mesh." << std::endl;
	exit ( 1 );                                                                        }

-------------------------------------------------------------------------------------------


void CellIterator::Over::TwoVerticesOfSeg::reset ( ) // virtual from CellIterator::Core
{	this->passage = 0;  }

void CellIterator::Over::TwoVerticesOfSeg::reset ( Cell::Core * cll )
// virtual from CellIterator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot reset an iterator this way." << std::endl;
	exit ( 1 );                                                      }

void CellIterator::Over::TwoVerticesOfSeg::advance ( )
// virtual from CellIterator::Core
#IFNDEF NDEBUG
{	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }  }
#ENDIF
{	this->passage++;  }

bool CellIterator::Over::TwoVerticesOfSeg::in_range ( )
// virtual from CellIterator::Core
{	return this->passage < 2;   }

Cell::Core * CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first base (negative) then tip (positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage == 0 )
	{	assert ( this->seg_p->base_p );
		return this->seg_p->base_p;     }  // negative vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }  // positive vertex
}

Cell::Core * CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first tip (positive) then base (negative)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage )
	{	assert ( this->seg_p->base_p );
		return this->seg_p->base_p;     }    // negative vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }  } // positive vertex


Cell::Core * CellIterator::Over::TwoVerticesOfPosSeg::NormalOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first base then tip (both positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage == 0 )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }     // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }             }  // positive vertex


Cell::Core * CellIterator::Over::TwoVerticesOfPosSeg::ReverseOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first tip then base (both positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }      // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }              }  // positive vertex


Cell::Core * CellIterator::Over::TwoVerticesOfNegSeg::NormalOrder::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first base (negative) then tip (positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }      // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		assert ( this->seg_p->tip_p->reverse_p );
		return this->seg_p->tip_p->reverse_p;      }   }  // negative vertex
	

Cell::Core * CellIterator::Over::TwoVerticesOfNegSeg::ReverseOrder::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first tip (positive) then base (negative)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage == 0 )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }      // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		assert ( this->seg_p->tip_p->reverse_p );
		return this->seg_p->tip_p->reverse_p;      }   }  // negative vertex


Cell::Core * CellIterator::Over::TwoVerticesOfNegSeg::NormalOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first base then tip (both positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }        // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }                }  // positive vertex


Cell::Core * CellIterator::Over::TwoVerticesOfNegSeg::ReverseOrder::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	// iterate over the two vertices, first tip then base (both positive)
#IFNDEF NDEBUG
	if ( this->passage >= 2 )
	{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
							<< __extension__ __PRETTY_FUNCTION__ << ": ";
		std::cout << "This iterator is out of range." << std::endl;
		exit ( 1 );                                                  }
#ENDIF
	if ( this->passage == 0 )
	{	assert ( this->seg_p->base_p );
		assert ( this->seg_p->base_p->reverse_p );
		return this->seg_p->base_p->reverse_p;     }      // positive vertex
	else
	{	assert ( this->seg_p->tip_p );
		return this->seg_p->tip_p;      }              }  // positive vertex


void CellIterator::Over::CellsOfFuzzyMesh::reset ( ) // virtual from CellIterator::Core
{	this->iter = this->list_p->begin();  }

void CellIterator::Over::CellsOfFuzzyMesh::reset ( Cell::Core * cll )
// virtual from CellIterator::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
						<< __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot reset an iterator this way." << std::endl;
	exit ( 1 );                                                      }

void CellIterator::Over::CellsOfFuzzyMesh::advance ( )
// virtual from CellIterator::Core
{	this->iter++;  }

bool CellIterator::Over::CellsOfFuzzyMesh::in_range ( )
// virtual from CellIterator::Core
{	return this->iter != this->list_p->end();  }

Cell::Core * CellIterator::Over::CellsOfFuzzyMesh::AsTheyAre::deref ( )
// virtual from CellIterator::Core
{	return * ( this->iter );  }

Cell::Core * CellIterator::Over::CellsOfFuzzyMesh::ReverseEachCell::deref ( )
// virtual from CellIterator::Core
{	Cell::Core * cll_p = ( * this->iter )->reverse_p;
	assert ( cll_p );  return cll_p;                   }

Cell::Core * CellIterator::Over::CellsOfFuzzyMesh::ForcePositive::deref ( )
// virtual from CellIterator::Core
{	return ( * this->iter )->get_positive();  }


	
