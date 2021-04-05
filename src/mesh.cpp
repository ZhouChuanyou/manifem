
// mesh.cpp 2021.04.05

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

#include <forward_list>

#include "mesh.h"
#include "iterator.h"

using namespace maniFEM;


size_t Mesh::maximum_dimension_plus_one { 4 };  // static data member

// we keep here the topological dimension of the largest mesh we intend to build
// '4' means three-dimensional meshes (cubes, tetrahedra, etc)
// '3' means two-dimensional meshes, including surfaces in R^3
// '2' would be for just polygonal lines
// '1' doesn't make much sense - just points ?

// see method Mesh::set_max_dim and paragraph 9.5 in the manual

// static data members :
std::vector < size_t > Cell::double_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0. );
std::vector < size_t > Cell::double_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0. );
std::vector < size_t > Cell::size_t_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::size_t_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::short_int_heap_size_pos ( Mesh::maximum_dimension_plus_one, 0 );
std::vector < size_t > Cell::short_int_heap_size_neg ( Mesh::maximum_dimension_plus_one, 0 );

	
Cell::Core * const Cell::ghost { new Cell::NegativeVertex ( tag::ghost ) };
// static data member, see paragraph 9.14 in the manual

//-----------------------------------------------------------------------------//


Mesh::Mesh ( const tag::DeepCopy &, const Mesh & msh )
:	core { nullptr }, meth { & Mesh::Positive::methods_pos }
{	size_t d = msh.dim();
	if ( d > 1 )
		this->core = new Mesh::Positive ( tag::of_dimension, d+1, tag::minus_one );
	else
	{	assert ( d == 1 );
		this->core = new Mesh::OneDim::Positive ( );  }
	std::list < Cell::Core* > & l = msh.core->cells[d];
	std::list<Cell::Core*>::iterator it = l.begin();
	for ( ; it != l.end(); it++ )
	{	Cell::Core * cll_p = *it;
		Cell cll ( tag::whose_core_is, cll_p );
		cll.add_to ( *this );   }                                                    }

//	CellIterator it = msh.iter_over ( tag::cells_of_dim, d );
//	for ( it.reset(); it.in_range(); it++ )
//	{	Cell cll = *it;
	
//-----------------------------------------------------------------------------//


bool Cell::Positive::is_positive ( ) const  // virtual from Cell::Core
{	return true;  }

bool Cell::Core::Negative::is_positive ( ) const  // virtual from Cell::Core
{	return false;  }

bool Mesh::return_true ( ) // static
{	return true;  }

bool Mesh::return_false ( ) // static
{	return false;  }

Cell::Positive * Cell::Positive::get_positive ( )  // virtual from Cell::Core
{	return this;  }

Cell::Positive * Cell::Core::Negative::get_positive ( )  // virtual from Cell::Core
{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	return (Cell::Positive*) this->reverse_p;     }


Cell::Core * Cell::Negative::build_reverse ( )
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "do not use build_reverse on a negative cell" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Cell::PositiveVertex::build_reverse ( )
// virtual from Cell::Core
{	this->reverse_p = new Cell::NegativeVertex ( tag::reverse_of, this );
	return this->reverse_p;                                                    }

Cell::Core * Cell::PositiveSegment::build_reverse ( )
// virtual from Cell::Core
{	this->reverse_p = new Cell::NegativeSegment ( tag::reverse_of, this );
	return this->reverse_p;                                                     }

Cell::Core * Cell::PositiveHighDim::build_reverse ( )
// virtual from Cell::Core
{	this->reverse_p = new Cell::Negative ( tag::reverse_of, this );
	return this->reverse_p;                                                   }

Cell::Core * Cell::Positive::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core
{	if ( this->reverse_p == nullptr )  this->reverse_p = this->build_reverse();
	assert ( this->reverse_p );
	assert ( not this->reverse_p->is_positive() );
	return this->reverse_p;                                                      }

Cell::Core * Cell::Negative::reverse ( const tag::BuildIfNotExists & )
// virtual from Cell::Core
{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	return this->reverse_p;                      }

//-----------------------------------------------------------------------------//


size_t Cell::PositiveVertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::NegativeVertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::PositiveSegment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::NegativeSegment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::PositiveHighDim::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->boundary_p );
	return this->boundary_p->get_dim_plus_one();  }

size_t Cell::NegativeHighDim::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->reverse_p );
	return this->reverse_p->get_dim();  }


size_t Mesh::ZeroDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 1;  }

size_t Mesh::Connected::OneDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 2;  }

size_t Mesh::Connected::HighDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->nb_of_cells.size();  }

size_t Mesh::MultiplyConnected::OneDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return 2;  }

size_t Mesh::MultiplyConnected::HighDim::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->nb_of_cells.size();  }

size_t Mesh::Fuzzy::get_dim_plus_one ( )  // virtual from Mesh::Core
{	return this->cells.size();  }

	
size_t Mesh::ZeroDim::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	return 2;  }
	
size_t Mesh::ZeroDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "zero-dimensional meshes have have no segments" << std::endl;
	exit ( 1 );                                                                                     }
	
size_t Mesh::ZeroDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d == 0 );  return 2;  }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::Vertices & )
// virtual from Mesh::Core
{	if ( first_ver == last_ver )  // closed loop
		return this->nb_of_segments;
	return this->nb_of_segments + 1;  }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	return this->nb_of_segments;  }
	
size_t Mesh::Connected::OneDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	if ( d == 1 ) return this->nb_of_segments;
	assert ( d == 0 );
	if ( first_ver == last_ver )  // closed loop
		return this->nb_of_segments;
	return this->nb_of_segments + 1;             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d < this->get_dim_plus_one() );
	assert ( this->nb_of_cells.size() > d );
	return this->nb_of_cells[d];             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::Vertices )
// virtual from Mesh::Core
{	assert ( this->nb_of_cells.size() > 0 );
	return this->nb_of_cells[0];             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	assert ( this->nb_of_cells.size() > 1 );
	return this->nb_of_cells[1];             }
	
size_t Mesh::Fuzzy::number_of ( const tag::CellsOfDim &, const size_t d )
// virtual from Mesh::Core
{	assert ( d < this->get_dim_plus_one() );
	assert ( this->nb_of_cells.size() > d );
	return this->nb_of_cells[d];             }
	
size_t Mesh::Fuzzy::number_of ( const tag::Vertices )
// virtual from Mesh::Core
{	assert ( this->cells.size() > 0 );
	return this->cells[0].size();       }
	
size_t Mesh::Fuzzy::number_of ( const tag::Segments & )
// virtual from Mesh::Core
{	assert ( 1 < this->get_dim_plus_one() );
	assert ( this->cells.size() > 1 );
	return this->cells[1].size();       }
	

Cell::Core * Mesh::Core::first_vertex ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::last_vertex ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last vertex" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::first_segment ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have first segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::last_segment ( )  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "only one-dimensional meshes have last segment" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Connected::OneDim::first_vertex ( )
// virtual from Mesh::Core, here overriden
// returns a negative vertex
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	return this->first_ver;                         }

Cell::Core * Mesh::Connected::OneDim::last_vertex ( )
// virtual from Mesh::Core, here overriden
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	return this->last_ver;                          }

Cell::Core * Mesh::Connected::OneDim::first_segment ( )
// virtual from Mesh::Core, here overriden
{	this->order();
	assert ( this->first_ver );
	assert ( this->first_ver != Cell::ghost );
	assert ( not this->first_ver->is_positive() );
	assert ( this->last_ver );
	assert ( this->last_ver->is_positive() );
	Cell::Core * neg_ver = this->first_ver;
	assert ( neg_ver );  assert ( not neg_ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = neg_ver->cell_behind_within.find ( this );
	assert ( it != neg_ver->cell_behind_within.end() );
	assert ( it->second );  // check what happens for an empty mesh
	return it->second;                                               }

Cell::Core * Mesh::Connected::OneDim::last_segment ( )
// virtual from Mesh::Core, here overriden
{	Cell::Core * ver = this->last_vertex();
	assert ( ver );  assert ( ver->is_positive() );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = ver->cell_behind_within.find ( this );
	assert ( it != ver->cell_behind_within.end() );
	assert ( it->second );  // check what happens for an empty mesh
	return it->second;                                               }

//-----------------------------------------------------------------------------//


Mesh Cell::PositiveVertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                                                 }

Mesh Cell::NegativeVertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                                                 }

Mesh Cell::PositiveSegment::boundary ( )  // virtual from Cell::Core
{	return Mesh ( tag::boundary_of, tag::positive, tag::segment, seg );  }

Mesh Cell::NegativeSegment::boundary ( )  // virtual from Cell::Core
{	assert ( this->reverse_p );
	Cell::PositiveSegment * rev_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( this->reverse_p );
	return Mesh ( tag::boundary_of, tag::positive, tag::segment, rev_seg,
                tag::reversed, tag::cells_surely_exist                  );  }

Mesh Cell::PositiveHighDim::boundary ( )  // virtual from Cell::Core
{	return Mesh ( tag::whose_core_is, this->boundary_p, tag::is_positive );   }

Mesh Cell::NegativeHighDim::boundary ( )  // virtual from Cell::Core
{	Cell::PositiveHighDim * cll = Util::assert_cast
		< Cell::Core*, Cell::PositiveHighDim* > ( this->reverse_p );
	// return reversed mesh, faces surely exist
	return Mesh ( tag::whose_core_is, cll->boundary_p,
	              tag::is_negative, tag::cells_surely_exist );        }   }

//-----------------------------------------------------------------------------//


// the two methods below are only relevant for STSI meshes
// so we forbid execution in Mesh::Core and then override them in Mesh::STSI

Cell::Core * Mesh::Core::cell_in_front_of  // virtual
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se                                                  ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_in_front as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }

Cell::Core * Mesh::Core::cell_behind  // virtual
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se                                                  ) const
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "cell_behind as seen_from : use on STSI meshes only" << std::endl;
	exit ( 1 );                                                                                     }


// the four methods below are only relevant for STSI meshes

Cell::Core * Mesh::STSI::cell_in_front_of
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::MayNotExist                                                        ) const

{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	Cell::Core * face_rev_p = face_p->reverse_p;
	assert ( face_rev_p );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = face_rev_p->cell_behind_within.find ( this->core );
	assert ( it != face_rev_p->cell_behind_within.end() );  //  to finish !
	Cell::Core * cll_rev_p = it->second;
	assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
	return cll_rev_p->reverse_p;                                        }  }

}	

Cell::Core * Mesh::STSI::cell_in_front_of
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se                                                  ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument

{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	Cell::Core * face_rev_p = face_p->reverse_p;
	assert ( face_rev_p );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = face_rev_p->cell_behind_within.find ( this->core );
	assert ( it != face_rev_p->cell_behind_within.end() );  //  to finish !
	Cell::Core * cll_rev_p = it->second;
	assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
	return cll_rev_p->reverse_p;                                        }  }

}	

Cell::Core * Mesh::STSI::cell_behind
// virtual from Mesh::Core, here overriden
( const Cell face, const tag::SeenFrom &, const Cell neighbour,
  const tag::MayNotExist                                        ) const
	
{	assert ( this->get_dim_plus_one() == face.dim() + 2 );
	std::map < Mesh::Core*, Cell::Core* >::const_iterator
		it = face.core->cell_behind_within.find ( this->core );
	if ( it == face.core->cell_behind_within.end() ) return nullptr;
		// nothing behind us, we are touching the boundary
	Cell::Core * cll_p = it->second;
	if ( cll_p ) return cll_p;
	// now cll_p is null which means face_p is a singular point,
	// i.e. a touching point or a self-intersection
	// so its neighbours are stored separately
	std::vector < std::pair < Cell, Cell > > ::iterator itt = this->singular.begin();
	for ( itt.reset(); itt.in_range(); itt++ );
	{	std::pair < Cell, Cell > & p = *itt;
		Cell c1 & = p.first;
		Cell c2 & = p.second;
		Cell face_rev = face.reverse ( tag::surely_exists );
		if ( face_rev.belongs_to ( c1.boundary(), tag::oriented ) ) return c2.core;
		if ( face_rev.belongs_to ( c2.boundary(), tag::oriented ) ) return c1.core;
		assert ( false );                                                            }    }


Cell::Core * Mesh::STSI::cell_behind
// virtual from Mesh::Core, here overriden
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se                                                  ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
	
{	assert ( this->get_dim_plus_one() == face_p->get_dim() + 2 );
	std::map<Mesh::Core*,Cell::Core*>::const_iterator
		it = face_p->cell_behind_within.find ( this->core );
	assert ( it != face_p->cell_behind_within.end() );
	Cell::Core * cll_p = it->second;
	if ( cll_p ) return cll_p;
	// now cll_p == nullptr which means face_p is a singular point,
	// i.e. a touching point or a self-intersection
	// so its neighbours are stored separately

}
	
//-----------------------------------------------------------------------------//


bool Cell::PositiveHighDim::belongs_to ( Mesh::Core * msh, const tag::Oriented & ) const
// virtual from Cell::Core

{	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap = this->meshes[0];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	if ( it == mmap.end() ) return false;
	const Cell::field_to_meshes & field = it->second;
	if ( field.counter_pos != 1 ) return false;
	assert ( field.counter_neg == 0 );
	return true;                                                                     }


bool Cell::PositiveHighDim::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		this->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                      }


bool Cell::Negative::belongs_to ( Mesh::Core * msh, const tag::Oriented & ) const
// virtual from Cell::Core

{	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	Cell::Positive * rev = static_cast < Cell::Positive * > ( this->reverse_p );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap = rev->meshes[0];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	if ( it == mmap.end() ) return false;
	const Cell::field_to_meshes & field = it->second;
	if ( field.counter_neg != 1 ) return false;
	assert ( field.counter_pos == 0 );
	return true;                                                                       }


bool Cell::Negative::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	Cell::Positive * rev = static_cast < Cell::Positive * > ( this->reverse_p );
	std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		rev->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                      }

//-----------------------------------------------------------------------------//


Cell::Core * Cell::Core::tip ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::Core::base ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::PositiveSegment::tip () { return this->tip_p;  }
// virtual from Cell::Core, overridden here

Cell::Core * Cell::PositiveSegment::base () { return this->base_p;  }
// virtual from Cell::Core, overridden here

Cell::Core * Cell::NegativeSegment::tip ()  // virtual from Cell::Core, overridden here
{	assert ( this->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( this->reverse_p );
	assert ( pos_seg->base_p );
	// assert ( pos_seg->base_p->reverse_p );
	return pos_seg->base_p->reverse_p;                                               }

Cell::Core * Cell::NegativeSegment::base ()  // virtual from Cell::Core, overridden here
{	assert ( this->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( this->reverse_p );
	assert ( pos_seg->tip_p );
	// assert ( pos_seg->tip_p->reverse_p );
	return pos_seg->tip_p->reverse_p;                                                }

//-----------------------------------------------------------------------------//

#ifndef NDEBUG


std::string Cell::PositiveHighDim::get_name ()  // virtual from Cell::Core
{	return this->name;  }

std::string Cell::Negative::get_name ()  // virtual from Cell::Core
{	return "r" + this->reverse_p->name;  }


void Cell::PositiveVertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is PositiveVertex " << this->name << std::endl;
	if ( this->meshes.size() > 0 )
	{	if ( this->meshes[0].size() > 0 )
				std::cout << "meshes of index 0, dim 0 (segments disguised as meshes)" << std::endl;
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
			it = this->meshes[0].begin(), it_e = this->meshes[0].end();
		for ( ; it != it_e; it++ )
		{	Cell::PositiveSegment * seg = (Cell::PositiveSegment*) it->first;
			std::cout << seg->name << " " << it->second.counter_pos
		                         << " " << it->second.counter_neg << "  ";  }
			if ( this->meshes[0].size() > 0 )  std::cout << std::endl;
		for ( size_t d = 1; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							itt = this->meshes[d].begin(), itt_e = this->meshes[d].end();
			for ( ; itt != itt_e; itt++ )
				std::cout << itt->first->get_name() << " " << itt->second.counter_pos
			                                      << " " << itt->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                        }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                             }


void Cell::NegativeVertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Cell::NegativeVertex, reverse of "
            << this->reverse_p->name << std::endl;            }


void Cell::PositiveSegment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is PositiveSegment " << this->name << std::endl;
	std::cout << "base : " << this->base_p->get_name() << std::endl;
	std::cout << "tip :  " << this->tip_p->get_name()  << std::endl;
	if ( this->meshes.size() > 0 )
	{	for ( size_t d = 0; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d << ", that is, dimension " << d+1 << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes[d].begin(), it_e = this->meshes[d].end();
			for ( ; it != it_e; it++ )
				std::cout << it->first->get_name() << " " << it->second.counter_pos
			                                     << " " << it->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                     }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                          }


void Cell::NegativeSegment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is NegativeSegment " << this->get_name() << std::endl;     }


void Cell::PositiveHighDim::print_everything ( )  // virtual from Cell::Core
	
{	size_t dim = this->get_dim();
	std::cout << "this is Cell::PositiveHighDim of dim "
	          << dim << " " << this->name << std::endl;
	if ( this->meshes.size() > 0 )
	{	std::cout << "meshes above me :" << std::endl;
		for ( size_t d = 0; d < this->meshes.size(); d++ )
		{	if ( this->meshes[d].size() > 0 )
				std::cout << "meshes of index " << d
			            << ", that is, dimension " << d+dim << std::endl;
			std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
							it = this->meshes[d].begin(), it_e = this->meshes[d].end();
			for ( ; it != it_e; it++ )
				std::cout << it->first->get_name() << " " << it->second.counter_pos
		                                       << " " << it->second.counter_neg << "  ";
			if ( this->meshes[d].size() > 0 )  std::cout << std::endl;                 }  }
	if ( this->reverse_p ) std::cout << "has reverse" << std::endl;                       }


void Cell::NegativeHighDim::print_everything ( )  // virtual from Cell::Core
	
{	std::cout << "this is Cell::NegativeHighDim "
						<< this->reverse_p->get_name() << std::endl;   }


std::string Mesh::Positive::get_name ( )  // virtual from Mesh::Core

{	if ( this->cell_enclosed ) return "bdry_of_" + this->cell_enclosed->name;
	else return this->name;                                                   }

std::string Mesh::OneDim::Positive::get_name ( )  // virtual from Mesh::Core

{	if ( this->cell_enclosed ) return "bdry_of_" + this->cell_enclosed->name;
	else return this->name;                                                   }

void Mesh::Positive::print_everything ( )  // virtual from Mesh::Core

{	std::cout << "this is Mesh::Positive " << this->get_name() << std::endl;
	for ( size_t d = 0; d < this->cells.size(); d++ )
	{	std::cout << "cells of dim " << d << " :" << std::endl;
		int counter = 0;
		std::list<Cell::Core*>::iterator it = this->cells[d].begin(),
		                                 it_e = this->cells[d].end();
		for ( ; it != it_e; it++, counter++ ) std::cout << (*it)->get_name() << " ";
		if ( counter > 0 )  std::cout << std::endl;                                   }  }

void Mesh::OneDim::Positive::print_everything ( )  // virtual from Mesh::Core

{	std::cout << "this is Mesh::OneDim::Positive " << this->get_name() << std::endl;
	for ( size_t d = 0; d < this->cells.size(); d++ )
	{	std::cout << "cells of dim " << d << " :" << std::endl;
		int counter = 0;
		std::list<Cell::Core*>::iterator it = this->cells[d].begin(),
		                                 it_e = this->cells[d].end();
		for ( ; it != it_e; it++, counter++ ) std::cout << (*it)->get_name() << " ";
		if ( counter > 0 )  std::cout << std::endl;                                   }  }

#endif

 //////////////////////////////////////////////////////////////////
////////////    add/remove a cell to/from a mesh    ///////////////
///////////////////////////////////////////////////////////////////
////////////     this is where the magic happens     //////////////
//////////// (and also where we get our hands dirty) //////////////
//////////////////////////////////////////////////////////////////


void Cell::PositiveVertex::glue_on_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::PositiveSegment::glue_on_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->add_to_seg ( this );
	// 'add_to_seg' is virtual, so the computer will choose the right version
	// (Cell::PositiveVertex::add_to_seg or Cell::NegativeVertex::add_to_seg)
	this->glue_common ( ver );       }
	

void Cell::PositiveHighDim::glue_on_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	face->add_to_bdry ( this->boundary_p );
	// 'add_to_bdry' is virtual, so the computer will choose the right version
	this->glue_common ( face );                         }


void Cell::Negative::glue_on_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_p );
	assert ( this->reverse_p );
	this->reverse_p->glue_on_my_bdry ( cll->reverse_p );  }


void Cell::PositiveVertex::cut_from_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::PositiveSegment::cut_from_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->remove_from_seg ( this );
	// 'remove_from_seg' is virtual, so the computer will choose the right version
	// (Cell::PositiveVertex::remove_from_seg or Cell::NegativeVertex::remove_from_seg)
	this->cut_common ( ver );       }
	

void Cell::PositiveHighDim::cut_from_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	face->remove_from ( this->boundary_p );
	// 'remove_from' is virtual, so the computer will choose the right version
	this->cut_common ( face );                            }


void Cell::Negative::cut_from_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_p );
	assert ( this->reverse_p );
	this->reverse_p->cut_from_my_bdry ( cll->reverse_p );  }

//-----------------------------------------------------------------------------//


namespace {  // anonymous namespace, mimiks static linkage

// here is where the low-level linking between cells and meshes happens
// ***_deep_connections and ***_deep_connections_rev
// are called from Mesh::***::add_*** and remove_***

	
inline void add_link_zero_dim  // hidden in anonymous namespace
( Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.

// this function updates ver->meshes_same_dim
// also, calls the virtual method msh->add_to_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes, returns garbage

{	assert ( ver );  assert ( seg );
///////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map < Cell::PositiveSegment *, short int > maptype;       //
	maptype & cmd = ver->segments;                                           //
	maptype::iterator lb = cmd.lower_bound(seg);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(seg,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(seg),                                        // 
	      std::forward_as_tuple(1,seg->add_to_cells(ver,0)) );               //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( ver->meshes_same_dim.find(msh) == ver->meshes_same_dim.end() );   //
//	msh->cells[ver->dim].push_front(o_ver);                                    //
//	Cell::field_to_meshes field;                                               //
//	field.counter_pos = cp;                                                    //
//	field.counter_neg = cn;                                                    //
//	field.where = msh->cells[ver->dim]->begin();                               //
//	ver->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_zero_dim


inline void add_link_zero_dim_rev  // hidden in anonymous namespace
( Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.

// this function updates ver->meshes_same_dim
// also, calls the virtual method msh->add_to_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes, returns garbage

{	assert ( ver );  assert ( seg );
///////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map < Cell::PositiveSegment *, short int > maptype;       //
	maptype & cmd = ver->segments;                                           //
	maptype::iterator lb = cmd.lower_bound(seg);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(seg,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(seg),                                        // 
	      std::forward_as_tuple(1,seg->add_to_cells(ver,0)) );               //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( ver->meshes_same_dim.find(msh) == ver->meshes_same_dim.end() );   //
//	msh->cells[ver->dim].push_front(o_ver);                                    //
//	Cell::field_to_meshes field;                                               //
//	field.counter_pos = cp;                                                    //
//	field.counter_neg = cn;                                                    //
//	field.where = msh->cells[ver->dim]->begin();                               //
//	ver->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_zero_dim_rev


inline void add_link_same_dim  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is positive.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 == cll_dim + 1 );
///////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;    //
	maptype & cmd = cll->meshes_same_dim;                                    //
	maptype::iterator lb = cmd.lower_bound(msh);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(msh),                                        // 
	      std::forward_as_tuple(1,msh->add_to_cells(cll,cll_dim)) );         //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(o_cll);                                    //
//	Cell::field_to_meshes field;                                               //
//	field.counter_pos = cp;                                                    //
//	field.counter_neg = cn;                                                    //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim


inline void add_link_same_dim_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll, Mesh::Core * msh )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is negative.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->add_to_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes, returns garbage

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( ! o_cll->is_positive() );
	assert ( ( sign == 1 ) or ( sign == -1 ) );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 == cll_dim + 1 );
///////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL         //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> maptype;    //
	maptype & cmd = cll->meshes_same_dim;                                    //
	maptype::iterator lb = cmd.lower_bound(msh);                             //
	assert ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) );     //
	cmd.emplace_hint ( lb, std::piecewise_construct,                         //
	      std::forward_as_tuple(msh),                                        // 
	      std::forward_as_tuple(-1,msh->add_to_cells(o_cll,cll_dim)) );      //
/////////  code below is conceptually equivalent to the above  //////////////////
//	assert ( cll->meshes_same_dim.find(msh) == cll->meshes_same_dim.end() );   //
//	msh->cells[cll->dim].push_front(o_cll);                                    //
//	Cell::field_to_meshes field;                                               //
//	field.counter_pos = cp;                                                    //
//	field.counter_neg = cn;                                                    //
//	field.where = msh->cells[cll->dim]->begin();                               //
//	cll->meshes_same_dim[msh] = field;                                         //
/////////////////////////////////////////////////////////////////////////////////

} // end of add_link_same_dim


inline void add_link  // hidden in anonymous namespace
( Cell::Positive * cll, Mesh::Core * msh, short int cp, short int cn )

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation.

// this function updates cll->meshes
// also, calls the virtual method msh->add_to_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// returning an iterator within msh->cells
// for other kinds of meshes, returns garbage

{	assert ( cll );  assert ( msh );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
/////////////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL           //
	typedef std::map<Mesh::Core*,Cell::field_to_meshes> maptype;               //
	maptype & cmd = cll->meshes[dif_dim];                                      //
	maptype::iterator lb = cmd.lower_bound(msh);                               //
	if ( ( lb == cmd.end() ) or ( cmd.key_comp()(msh,lb->first) ) )            //
		cmd.emplace_hint ( lb, std::piecewise_construct,                         //
		      std::forward_as_tuple(msh),                                        // 
		      std::forward_as_tuple(cp,cn,msh->add_to_cells(cll,cll_dim)) );     //
	else                                                                       //
	{ lb->second.counter_pos += cp;                                            //
	  lb->second.counter_neg += cn;   }                                        //
////////// code below is conceptually equivalent to the above /////////////////
//	if ( cll->meshes[dif_dim].find(msh)==cll->meshes[dif_dim].end() )   //
//	{ msh->cells[cll->dim].push_front(cll);                             //
//	  Cell::field_to_meshes field;                                      //
//	  field.counter_pos = cp;                                           //
//	  field.counter_neg = cn;                                           //
//	  field.where = msh->cells[cll->dim]->begin();                      //
//	  cll->meshes[dif_dim][msh] = field;             }                  //
//	else                                                                //
//	{ cll->meshes[dif_dim][msh].counter_pos += cp;                      //
//	  cll->meshes[dif_dim][msh].counter_neg += cn;   }                  //
//////////////////////////////////////////////////////////////////////////

} // end of add_link


inline void link_face_to_msh  // hidden in anonymous namespace
( Cell::Core * face, Mesh::Core * msh, short int cp, short int cn )
// just a block of code for make_deep_connections

{	assert ( face.is_positive() );
	assert ( face->get_dim() + 1 < msh->get_dim_plus_one() ):
	Cell::Positive * face_p = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	size_t dif_dim = msh->get_dim_plus_one() - face->get_dim() - 1;
	// implement a separate function link_cll_to_msh !!
	assert ( dif_dim < face_p->meshes.size() );
	std::map < Mesh::Core*, Cell::field_to_meshes > & cemd = face_p->meshes[dif_dim];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		map_iter = cemd.find(cell_bdry);
	assert ( map_iter != cemd.end() );
	Cell::field_to_meshes & mis = map_iter->second;
	add_link ( face_p, msh, cp*mis.counter_pos + cn*mis.counter_neg,
                          cn*mis.counter_pos + cp*mis.counter_neg );                }
	

inline void link_face_to_higher  // hidden in anonymous namespace
( Cell::Core * face, Cell::Positive * pmce, short int cp, short int cn )
// just a block of code for make_deep_connections

{	assert ( face );
	assert ( face->get_dim() + 2 < pmce->get_dim() ):
	Cell::Positive * face_p = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	Cell::Positive * pmcenv = Util::assert_cast
		< Cell::Core*, Cell::PositiveNotVertex* > ( pmce );
	// pmce->meshes[0] is empty, we use pmcenv->meshes_same_dim
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> & cemd = pmcenv->meshes_same_dim;
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim>::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() ==
						 dif_dim + this->cell_enclosed->get_dim() + 1 );
		size_t current_mesh_dim = map_iter->first->get_dim_plus_one() - 1;
		assert ( current_mesh_dim == mce->dim() + dif_dim );
		Cell::field_to_meshes & mis = map_iter->second;
		if ( mis.sign == 1 ) add_link ( face_p, map_iter.first, cp, cn );
		else
		{	assert ( mis.sign == -1 );
			add_link ( face_p, map_iter.first, cn, cp );  }                  }
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	std::map<Mesh::Core*,Cell::field_to_meshes> & cemd = pmce->meshes[dif_dim];
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() ==
							 dif_dim + this->cell_enclosed->get_dim() + 1 );
			size_t current_mesh_dim = map_iter->first->get_dim_plus_one() - 1;
			assert ( current_mesh_dim == mce->dim() + dif_dim );
			Cell::field_to_meshes & mis = map_iter->second;
			add_link ( face_p, map_iter.first,
		             cp*mis.counter_pos + cn*mis.counter_neg,
		             cn*mis.counter_pos + cp*mis.counter_neg );              }        }          }

//	Mesh::add_link ( lower_cell_tri.obj, lower_cell_tri.obj, msh
//	         cell_counter_pos*mesh_counter_pos + cell_counter_neg*mesh_counter_neg,
//	         cell_counter_pos*mesh_counter_neg + cell_counter_neg*mesh_counter_pos  );


inline void compute_sign  // hidden in anonymous namespace
( short int & cp, short int & cn, Cell face, Mesh::Core * cell_bdry )
// just a block of code for make_deep_connections

// 'face' is a face of 'cll', so it is a cell in 'cell_bdry' of maximum dimension
// we just want to know if it is positively oriented or not

{	assert ( face.is_positive() );
	assert ( face.dim() + 1 == cell_bdry->get_dim_plus_one() ):
	Cell::PositiveNotVertex * face_p = Util::assert_cast
		< Cell::Core*, Cell::PositiveNotVertex* > ( face.core );
	std::map < Mesh::Core*, Cell::field_to_meshes_same_dim > &
		cemd = face_p->meshes_same_dim;
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim>::iterator
		map_iter = cemd.find(cell_bdry);
	assert ( map_iter != cemd.end() );
	short int sign = map_iter->second.sign;
	if ( sign == 1 )  {  cp = 1;  cn = 0;  }
	else  {  assert ( sign == -1 );  cp = 0;  cn = 1;               }

	
inline void compute_cp_cn  // hidden in anonymous namespace
( short int & cp, short int & cn, Cell face, Mesh::Core * cell_bdry )
// just a block of code for make_deep_connections

// 'face' is a vertex or side of 'cll', so it is a cell in 'cell_bdry', of lower dimension
// we want to know how many times it appears as positive or as negative within 'cell_bdry'

{	assert ( face );
	assert ( face->get_dim() + 1 < cell_bdry->get_dim_plus_one() ):
	size_t dif_dim = cell_bdry->get_dim_plus_one() - face->get_dim() - 1;
	Cell::Positive * face_p = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	std::map < Mesh::Core*, Cell::field_to_meshes > & cemd = face_p->meshes[dif_dim];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		map_iter = cemd.find(cell_bdry);
	assert ( map_iter != cemd.end() );
	Cell::field_to_meshes & mis = map_iter->second;
	cp = mis.cp;  cn = mis.cn;                                                         }

	
inline void make_deep_connections_0d  // hidden in anonymous namespace
( Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// make far connections when adding a positive vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim ( ver, seg );
	link_face_to_higher ( ver, seg, 1, 0 );  }
	

inline void make_deep_connections_0d_rev  // hidden in anonymous namespace
( Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// make far connections when adding a negative vertex

{	assert ( ver );  assert ( seg );
	add_link_zero_dim_rev ( ver, seg );
	link_face_to_higher ( ver, seg, 0, 1 );  }
	

inline void make_deep_connections  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh, const tag::MeshIsNotBdry & )

// make far connections when adding a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case cll.get_dim() == 0
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( cll );
	assert ( msh->cell_enclosed == nullptr );

	add_link_same_dim ( cll, msh );
	
	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
	  link_face_to_msh ( face, msh, cp, cn );                                   }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			link_face_to_msh ( fface, msh, ccp, ccn );      }               }

} // end of make_deep_connections with tag::mesh_is_not_bdry


inline void make_deep_connections  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh, const tag::MeshIsBdry & )

// make far connections when adding a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case cll.get_dim() == 0
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim ( cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive * pmce = Util::assert_cast < Cell::Core*, Cell::Positive* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 1, 0 );

	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
	  link_face_to_msh ( face, msh, cp, cn );
	  link_face_to_higher ( face, pmce, cp, cn );                                    }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			link_face_to_msh ( fface, msh, ccp, ccn );
			link_face_to_higher ( fface, pmce, ccp, ccn );   }            }

} // end of make_deep_connections with tag::mesh_is_not_bdry


inline void make_deep_connections_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll,
  Mesh::Core * msh, const tag::MeshIsNotBdry &          )

// make far connections when adding a negative cell
// see paragraph 10.1 in the manual
	
{	assert ( cll != o_cll );
	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case cll.get_dim() == 0
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( msh->cell_enclosed == nullptr );

	add_link_same_dim_rev ( o_cll, cll, msh );
	
	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh'
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
		// we switch the two counters
	  link_face_to_msh ( face, msh, cn, cp );                                        }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh'
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			// we switch the two counters
			link_face_to_msh ( fface, msh, ccn, ccp );      }          }

} // end of make_deep_connections_rev with tag::mesh_is_not_bdry


inline void make_deep_connections_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll,
  Mesh::Core * msh, const tag::MeshIsBdry &             )

// make far connections when adding a negative cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// make_deep_connections_0d deals with the case cll.get_dim() == 0
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	add_link_same_dim_rev ( o_cll, cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive * pmce = Util::assert_cast < Cell::Core*, Cell::Positive* > ( mce );
	// link 'cll' to all meshes above 'this->cell_enclosed' (of all dimensions)
	link_face_to_higher ( cll, pmce, 0, 1 );

	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // add link from face to 'msh' and
		// to all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
		// we switch the two counters
	  link_face_to_msh ( face, msh, cn, cp );
	  link_face_to_higher ( face, pmce, cn, cp );                                    }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // add link from face to 'msh' and
			// to all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			// we switch the two counters
			link_face_to_msh ( fface, msh, ccn, ccp );
			link_face_to_higher ( fface, pmce, ccn, ccp );   }           }

} // end of make_deep_connections_rev with tag::mesh_is_bdry


inline void remove_link_same_dim  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh )

// This function "removes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is positive.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdm = cmd.find(msh);
	assert ( cmdm != cmd.end() );
	assert ( cmdm->second.sign >= 1 );
	msh->remove_from_cells ( cll, cll_dim, cmdm->second.where );
	cmd.erase(cmdm);

} // end of remove_link_same_dim


inline void remove_link_same_dim_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll, Mesh::Core * msh )

// This function "removes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here the dimensions are equal and the cell is negative.

// this function updates cll->meshes_same_dim
// also, calls the virtual method msh->remove_from_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 == cll_dim + 1 );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes_same_dim> maptype;
	maptype & cmd = cll->meshes_same_dim;
	maptype::iterator cmdm = cmd.find(msh);
	assert ( cmdm != cmd.end() );
	assert ( cmdm->second.sign >= -1 );
	msh->remove_from_cells ( o_cll, cll_dim, cmdm->second.where );
	cmd.erase(cmdm);

} // end of remove_link_same_dim


inline void remove_link  // hidden in anonymous namespace
( Cell::Positive * cll, Mesh::Core * msh, short int cp, short int cn )

// This function "removes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// Here 'cll' has lower dimension, so it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with positive orientation,
// while 'cn' says the same thing for negative orientation.

// this function updates cll->meshes
// also, calls the virtual method msh->remove_from_cells
// which, for fuzzy and stsi meshes, updates msh->cells
// removing cll with the aid of the iterator provided
// for other kinds of meshes, does nothing

{	assert ( cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim + 1 );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes> maptype;
	maptype & cmd = cll->meshes[dif_dim];
	maptype::iterator cmdm = cmd.find(msh);
	assert ( cmdm != cmd.end() );
	short int c_p = cmdm->second.counter_pos -= cp;
	short int c_n = cmdm->second.counter_neg -= cn;
	assert ( ( c_p >= 0 ) and ( c_n >= 0 ) );
	if ( ( c_p == 0 ) and ( c_n == 0 ) )
	{	msh->remove_from_cells(cll,cll_dim,cmdm->second.where);
		cmd.erase(cmdm);                                                 }

} // end of remove_link


inline void unlink_face_from_msh  // hidden in anonymous namespace
( Cell::Core * face, Mesh::Core * msh, short int cp, short int cn )
// just a block of code for make_deep_connections

{	assert ( face->get_dim() + 1 < msh->get_dim_plus_one() ):
	Cell::Positive * face_p = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	size_t dif_dim = msh->get_dim_plus_one() - face->get_dim() - 1;
	assert ( dif_dim < face_p->meshes.size() );
	// implement a separate function unlink_link_cll_from_msh !!
	std::map < Mesh::Core*, Cell::field_to_meshes > & cemd = face_p->meshes[dif_dim];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		map_iter = cemd.find(cell_bdry);
	assert ( map_iter != cemd.end() );
	Cell::field_to_meshes & mis = map_iter->second;
	Mesh::remove_link ( face_p, msh, cp*mis.counter_pos + cn*mis.counter_neg,
                                   cn*mis.counter_pos +  cp*mis.counter_neg );        }
	

inline void unlink_face_from_higher
( Cell::Core * face, Cell::Positive * pmce, short int cp, short int cn )
// hidden in anonymous namespace
// just a block of code for make_deep_connections

{	assert ( face->get_dim() + 2 < pmce->get_dim() ):
	Cell::Positive * face_p = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( face );
	Cell::Positive * pmcenv = Util::assert_cast
		< Cell::Core*, Cell::PositiveNotVertex* > ( pmce );
	// pmce->meshes[0] is empty, we use pmcenv->meshes_same_dim
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim> & cemd = pmcenv->meshes_same_dim;
	std::map<Mesh::Core*,Cell::field_to_meshes_same_dim>::iterator map_iter;
	// we now loop over all meshes of given dimension
	for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
	{	assert ( map_iter->first->get_dim_plus_one() ==
						 dif_dim + this->cell_enclosed->get_dim() + 1 );
		size_t current_mesh_dim = map_iter->first->get_dim_plus_one() - 1;
		assert ( current_mesh_dim == mce->dim() + dif_dim );
		Cell::field_to_meshes & mis = map_iter->second;
		if ( mis.sign == 1 )
			Mesh::remove_link ( face_p, map_iter.first, cp, cn );
		else
		{	assert ( mis.sign == -1 );
			Mesh::remove_link ( face_p, map_iter.first, cn, cp );  }            }
	for ( size_t dif_dim = 1; dif_dim < pmce->meshes.size(); dif_dim++ )
	{	std::map<Mesh::Core*,Cell::field_to_meshes> & cemd = pmce->meshes[dif_dim];
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator map_iter;
		// we now loop over all meshes of given dimension
		for ( map_iter = cemd.begin(); map_iter != cemd.end(); ++map_iter )
		{	assert ( map_iter->first->get_dim_plus_one() ==
							 dif_dim + this->cell_enclosed->get_dim() + 1 );
			size_t current_mesh_dim = map_iter->first->get_dim_plus_one() - 1;
			assert ( current_mesh_dim == mce->dim() + dif_dim );
			Cell::field_to_meshes & mis = map_iter->second;
			Mesh::remove_link ( face_p, map_iter.first,
		                      cp*mis.counter_pos + cn*mis.counter_neg,
		                      cn*mis.counter_pos + cp*mis.counter_neg );      }  }  }


inline void break_deep_connections_0d  // hidden in anonymous namespace
( Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// break far connections when adding a positive vertex

{	assert ( ver );  assert ( seg );
	assert ( o_ver->is_positive() );
	remove_link_zero_dim ( ver, seg );
	unlink_face_from_higher ( ver, seg, 1, 0 );  }
	

inline void break_deep_connections_0d_rev  // hidden in anonymous namespace
( Cell::Core * o_ver, Cell::PositiveVertex * ver, Cell::PositiveSegment * seg )

// break far connections when adding a negative vertex

{	assert ( ver );  assert ( seg );
	assert ( o_ver );  assert ( ver != o_ver );
	remove_link_zero_dim_rev ( ver, seg );
	unlink_face_from_higher ( ver, seg, 0, 1 );   }
	

inline void break_deep_connections  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh, const tag::MeshIsNotBdry & )

// break far connections when removing a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 1 );
	// the case cll.get_dim() == 0 is dealt with separately,
	// see Mesh::ZeroDim::break_deep_connections
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( msh->cell_enclosed == nullptr );

	remove_link_same_dim ( cll, msh );
	
	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // unlink 'face' from 'msh'
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
	  unlink_face_from_msh ( face, msh, cp, cn );                                        }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // unlink 'face' from 'msh'
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			unlink_face_from_msh ( fface, msh, ccp, ccn );      }                               }

} // end of break_deep_connections with tag::mesh_is_not_bdry


inline void break_deep_connections  // hidden in anonymous namespace
( Cell::PositiveNotVertex * cll, Mesh::Core * msh, const tag::MeshIsBdry & )

// break far connections when removing a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( msh );
	assert ( msh->get_dim_plus_one() > 1 );
	// the case cll.get_dim() == 0 is dealt with separately,
	// see Mesh::ZeroDim::break_deep_connections
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	remove_link_same_dim ( cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive * pmce = Util::assert_cast < Cell::Core*, Cell::Positive* > ( mce );
	// unlink 'cll' from all meshes above 'this->cell_enclosed' (of all dimensions)
	unlink_face_from_higher ( cll, pmce, 1, 0 );

	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // unlink 'face' from 'msh' and
		// from all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
	  unlink_face_from_msh ( face, msh, cp, cn );
	  unlink_face_from_higher ( face, pmce, cp, cn );                                    }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // unlink 'face' from 'msh' and
			// from all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			link_face_to_msh ( fface, msh, ccp, ccn );
			link_face_to_higher ( fface, pmce, ccp, ccn );   }                              }

} // end of break_deep_connections with tag::mesh_is_bdry


inline void break_deep_connections_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll,
  Mesh::Core * msh, const tag::MeshIsNotBdry &          )

// break far connections when removing a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// the case cll.get_dim() == 0 is dealt with separately,
	// see Mesh::ZeroDim::break_deep_connections_rev
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );
	assert ( msh->cell_enclosed == nullptr );

	remove_link_same_dim_rev ( cll, msh );
	
	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // unlink 'face' from 'msh'
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
		// we switch the two counters
	  unlink_face_from_msh ( face, msh, cn, cp );                                        }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell::Core * fface = (*itt).core;  // unlink 'face' from 'msh'
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			// we switch the two counters
			unlink_face_from_msh ( fface, msh, ccn, ccp );      }        }

} // end of break_deep_connections_rev with tag::mesh_is_not_bdry


inline void break_deep_connections_rev  // hidden in anonymous namespace
( Cell::Core * o_cll, Cell::PositiveNotVertex * cll,
  Mesh::Core * msh, const tag::MeshIsBdry &            )

// break far connections when removing a positive cell
// see paragraph 10.1 in the manual
	
{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll != o_cll );
	assert ( msh->get_dim_plus_one() > 1 );
	// the case cll.get_dim() == 0 is dealt with separately,
	// see Mesh::ZeroDim::break_deep_connections_rev
	size_t cll_dim = cll->get_dim();
	size_t cll_dim_m1 = Util::assert_diff ( cll_dim, 1 );
	assert ( msh->get_dim_plus_one() == cll_dim + 1 );

	remove_link_same_dim_rev ( cll, msh );
	
	// for all meshes strictly above msh
	Cell::Core * mce = msh->cell_enclosed;
	assert ( mce );
	Cell::Positive * pmce = Util::assert_cast < Cell::Core*, Cell::Positive* > ( mce );
	// unlink 'cll' from all meshes above 'this->cell_enclosed' (of all dimensions)
	unlink_cll_from_higher ( cll, pmce, 0, 1 );

	Mesh cll_bdry = cll->boundary();
	CellIterator it = cll_bdry.iterator
		( tag::over_cells, tag::of_max_dim, d, tag::force_positive );
	// talvez implementar um iterador especial que devolva cp e cn
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face = (*it).core;  // unlink 'face' from 'msh' and
		// from all meshes above 'msh->cell_enclosed' (of all dimensions)
		short int cp, cn;
		compute_sign ( cp, cn, face, cell_bdry );
		assert ( ( ( cp == 1 ) and ( cp == 0 ) ) or ( ( cp == 0 ) and ( cp == 1 ) ) );
		// we switch the two counters
	  unlink_face_from_msh ( face, msh, cn, cp );
	  unlink_face_from_higher ( face, pmce, cn, cp );                                    }
	for ( size_t d = 0; d < cll_dim_m1 ; d++ )
	{	CellIterator itt = cll_bdry.iterator  // as they are : positive
			( tag::over_cells, tag::of_dim, d, tag::as_they_are );
		// talvez implementar um iterador especial que devolva cp e cn
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell:Core * fface = (*itt).core;  // unlink 'face' from 'msh' and
			// from all meshes above 'msh->cell_enclosed' (of all dimensions)
			short int ccp, ccn;
			compute_cp_cn ( ccp, ccn, fface, cell_bdry );
			// we switch the two counters
			link_face_to_msh ( fface, msh, ccn, ccp );
			link_face_to_higher ( fface, pmce, ccn, ccp );   }            }

} // end of break_deep_connections_rev with tag::mesh_is_bdry
	

void Cell::Core::add_to_seg ( Cell::PositiveSegment * seg )
// virtual, overriden by Cell::***tive::Vertex
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Core::remove_from_seg ( Cell::PositiveSegment * seg )
// virtual, overriden by Cell::***tive::Vertex
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }


void Cell::PositiveVertex::add_to_seg ( Cell::PositiveSegment * seg )
// virtual from Cell::Core, here overriden
	
{	assert ( seg );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	std::map < Cell::PositiveSegment *, short int > & tm0 = this->segments;
	// this->segments contains segments instead of zero-dimensional meshes
	assert ( tm0.find(seg) == tm0.end() );
	assert ( seg->tip_p == nullptr );
	seg->tip_p = this;
	tm0.emplace ( std::piecewise_construct,
		std::forward_as_tuple(seg), std::forward_as_tuple(1) );
	make_deep_connections_0d ( this, seg );                                     }


void Cell::PositiveVertex::remove_from_seg ( Cell::PositiveSegment * seg )
// virtual from Cell::Core, here overriden
	
{	assert ( seg );
	// assert that 'this' vertex belongs to the boundary of 'seg'
	std::map < Cell::PositiveSegment *, short int > & tm0 = this->segments;
	// this->segments contains segments instead of zero-dimensional meshes
	std::map<Cell::PositiveSegment*,short int>::iterator it = tm0.find ( seg );
	assert ( it != tm0.end() );
	assert ( seg->tip_p == this );
	seg->tip_p = nullptr;
	assert ( it->second == 1 );
	tm0.erase ( it );
	break_deep_connections_0d ( this, seg );                                       }


void Cell::NegativeVertex::add_to_seg ( Cell::PositiveSegment * seg )
// virtual from Cell::Core, here overriden

{ assert ( seg );
	assert ( this->reverse_p );
	Cell::PositiveVertex * pos_ver = Util::assert_cast
		< Cell::Core*, Cell::PositiveVertex* > ( this->reverse_p );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	std::map < Cell::PositiveSegment *, short int > & pvm0 = pos_ver->segments;
	// this->segments contains segments instead of zero-dimensional meshes
	assert ( pvm0.find(seg) == pvm0.end() );
	assert ( seg->base_p == nullptr );
	seg->base_p = this;
	pvm0.emplace ( std::piecewise_construct,
		std::forward_as_tuple(seg), std::forward_as_tuple(-1) );
	make_deep_connections_0d_rev ( this, pos_ver, seg );                          }


void Cell::NegativeVertex::remove_from_seg ( Cell::PositiveSegment * seg )
// virtual from Cell::Core, here overriden

{ assert ( seg );
	assert ( this->reverse_p );
	Cell::PositiveVertex * pos_ver = Util::assert_cast
		< Cell::Core*, Cell::PositiveVertex* > ( this->reverse_p );
	// assert that 'this' vertex belongs to the mesh 'msh'
	std::map < Cell::PositiveSegment *, short int > & pvm0 = pos_ver->segments;
	// this->segments contains segments instead of zero-dimensional meshes
	std::map<Cell::PositiveSegment*,short int>::iterator it = pvm0.find ( seg );
	assert ( it != pvm0.end() );
	assert ( seg->base_p );
	seg->base_p = nullptr;
	assert ( it->second == -1 );
	pvm0.erase ( it );
	break_deep_connections_0d_rev ( this, pos_ver, seg );                           }


void Cell::PositiveVertex::add_to_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                  }


void Cell::PositiveVertex::remove_from_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::NegativeVertex::add_to_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                 }


void Cell::NegativeVertex::remove_from_mesh ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::PositiveSegment::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::PositiveSegment::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::NegativeSegment::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::NegativeSegment::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_seg ( this, tag::mesh_is_not_bdry );  }


void Cell::PositiveHighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::PositiveHighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::NegativeHighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::NegativeHighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_hd_cell ( this, tag::mesh_is_not_bdry );  }


void Cell::PositiveVertex::add_to_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::PositiveVertex::remove_from_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::NegativeVertex::add_to_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use add_to_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::NegativeVertex::remove_from_bdry ( Mesh::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Please use remove_from_seg." << std::endl;
	exit ( 1 );                                               }


void Cell::PositiveSegment::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_seg ( this, tag::mesh_is_bdry );  }


void Cell::PositiveSegment::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_seg ( this, tag::mesh_is_bdry );  }


void Cell::NegativeSegment::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_seg ( this, tag::mesh_is_bdry );  }


void Cell::NegativeSegment::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_seg ( this, tag::mesh_is_bdry );  }


void Cell::PositiveHighDim::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::PositiveHighDim::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::NegativeHighDim::add_to_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_hd_cell ( this, tag::mesh_is_bdry );  }


void Cell::NegativeHighDim::remove_from_bdry ( Mesh::Core * msh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_hd_cell ( this, tag::mesh_is_bdry );  }

//-----------------------------------------------------------------------------//


std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
Mesh::Core::add_to_cells ( Cell::Core * cll, const size_t d )
// virtual, here returns garbage, overriden by Mesh::Fuzzy and later by Mesh::STSI

// called from add_link_same_dim and add_link (both hidden in anonymous namespace below)
	
{	return static_cast
		< std::map<Mesh::Core*,Cell::field_to_meshes>::iterator > ( nullptr );  }


std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
Mesh::Fuzzy::add_to_cells ( Cell::Core * cll, const size_t d )
// virtual from Cell::Core, overriden here, later overriden again by Mesh::STSI

// called from add_link_same_dim and add_link (both hidden in anonymous namespace below)

// add a cell to 'this->cells[d]' list, return iterator into that list
		
{	assert ( d == cll->dim() );
	assert ( d < this->get_dim_plus_one() );
	std::list <Cell::Core*> & mcd = this->cells[d];
	mcd.push_front ( cll );
	return mcd.begin();                              }


std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
Mesh::STSI::add_to_cells ( Cell::Core * cll, const size_t d )
// virtual from Cell::Core, overriden here a second time

// called from add_link_same_dim and add_link (both hidden in anonymous namespace below)

// add a cell to 'this->cells[d]' list, return iterator into that list
	
{	return Mesh::Fuzzy::add_to_cells ( cll, d );  }


void Mesh::Core::remove_from_cells
( Cell::Core *, const size_t, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator )
// virtual, here does nothing, overriden by Mesh::Fuzzy and later by Mesh::STSI	
// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace below)
{	}


void Mesh::Fuzzy::remove_from_cells
( Cell::Core * cll, const size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator it )
// virtual from Cell::Core, overriden here, later overriden again by Mesh::STSI

// called from remove_link_same_dim and remove_link (both hidden in anonymous namespace below)

// remove a cell from 'this->cells[d]' list using the provided iterator

{	assert ( d == cll->dim() );
	assert ( d < this->get_dim_plus_one() );
	assert ( w != this->cells[d].end() );
	this->cells[d].erase(w);                   }


void Mesh::STSI::remove_from_cells
( Cell::Core * cll, const size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator it )
// virtual from Cell::Core, overriden here a second time

// remove a cell from 'this->cells[d]' list using the provided iterator

{	Mesh::Fuzzy::remove_from_cells ( cll, d, it );  }

//-----------------------------------------------------------------------------//


// see paragraph 10.1 in the manual

void Mesh::ZeroDim::add_pos_seg ( Cell::PositiveSegment *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_pos_seg ( Cell::PositiveSegment *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_pos_seg ( Cell::PositiveSegment *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_pos_seg ( Cell::PositiveSegment *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_neg_seg ( Cell::NegativeSegment *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_neg_seg ( Cell::NegativeSegment *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_neg_seg ( Cell::NegativeSegment *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_neg_seg ( Cell::NegativeSegment *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_pos_hd_cell ( Cell::PositiveHighDim *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::add_pos_hd_cell ( Cell::PositiveHighDim *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::remove_pos_hd_cell ( Cell::PositiveHighDim *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::remove_pos_hd_cell ( Cell::PositiveHighDim *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::add_neg_hd_cell ( Cell::NegativeHighDim *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::add_neg_hd_cell ( Cell::NegativeHighDim *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::remove_neg_hd_cell ( Cell::NegativeHighDim *, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }

void Mesh::ZeroDim::remove_neg_hd_cell ( Cell::NegativeHighDim *, const tag::MeshIsBdry & )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }


void Mesh::Fuzzy::add_pos_seg ( Cell::PositiveSegment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) == seg->meshes_same_dim.end() );

	make_deep_connections ( seg, this, tag::mesh_is_not_bdry );

	assert ( seg->base_p->cell_behind_within.find(this) ==
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within[msh] = seg;
	assert ( seg->tip_p->cell_behind_within.find(this) ==
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within[this] = seg;                                 }


void Mesh::Fuzzy::add_pos_seg ( Cell::PositiveSegment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) == seg->meshes_same_dim.end() );

	make_deep_connections ( seg, this, tag::mesh_is_bdry );

	assert ( seg->base_p->cell_behind_within.find(this) ==
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within[msh] = seg;
	assert ( seg->tip_p->cell_behind_within.find(this) ==
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within[this] = seg;                                 }


void Mesh::Fuzzy::remove_pos_seg ( Cell::PositiveSegment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) != seg->meshes_same_dim.end() );

	break_deep_connections ( seg, this, tag::mesh_is_not_bdry );

	assert ( seg->base_p->cell_behind_within.find(this) !=
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within.erase(this);
	assert ( seg->tip_p->cell_behind_within.find(this) !=
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within.erase(this);                                 }
	

void Mesh::Fuzzy::remove_pos_seg ( Cell::PositiveSegment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) != seg->meshes_same_dim.end() );

	break_deep_connections ( seg, this, tag::mesh_is_bdry );

	assert ( seg->base_p->cell_behind_within.find(this) !=
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within.erase(this);
	assert ( seg->tip_p->cell_behind_within.find(this) !=
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within.erase(this);                                 }
	

void Mesh::Fuzzy::add_neg_seg ( Cell::NegativeSegment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( seg->reverse_p );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( seg_pos->meshes_same_dim.find(this) == seg_pos->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	make_deep_connections_rev ( seg, pos_seg, this, tag::mesh_is_not_bdry );

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	pos_seg->base_p->reverse_p->cell_behind_within[this] = seg;
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	pos_seg->tip_p->reverse_p->cell_behind_within[this] = seg;                            }

	
void Mesh::Fuzzy::add_neg_seg ( Cell::NegativeSegment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( seg->reverse_p );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( seg_pos->meshes_same_dim.find(this) == seg_pos->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	make_deep_connections_rev ( seg, pos_seg, this, tag::mesh_is_bdry );

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	pos_seg->base_p->reverse_p->cell_behind_within[this] = seg;
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	pos_seg->tip_p->reverse_p->cell_behind_within[this] = seg;                            }

	
void Mesh::Fuzzy::remove_neg_seg ( Cell::NegativeSegment * seg, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( seg->reverse_p );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim.find(this) != pos_seg->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	break_deep_connections_rev ( seg, pos_seg, this, tag::mesh_is_not_bdry );

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->base_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->base_p->reverse_p->cell_behind_within.erase(this);
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->tip_p->reverse_p->cell_behind_within.erase(this);                           }
	
	
void Mesh::Fuzzy::remove_neg_seg ( Cell::NegativeSegment * seg, const tag::MeshIsBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::PositiveSegment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::PositiveSegment* > ( seg->reverse_p );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim.find(this) != pos_seg->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	break_deep_connections_rev ( seg, pos_seg, this, tag::mesh_is_bdry );

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->base_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->base_p->reverse_p->cell_behind_within.erase(this);
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->tip_p->reverse_p->cell_behind_within.erase(this);                           }
	
	
void Mesh::Fuzzy::add_pos_hd_cell ( Cell::PositiveHighDim * cll, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) == cll->meshes_same_dim.end() );
	
	make_deep_connections ( cll, this, tag::mesh_is_not_bdry );
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		assert ( face_p->cell_behind_within.find(this) ==
		         face_p->cell_behind_within.end()         );
		face_p->cell_behind_within[this] = cll;              }                          }


void Mesh::Fuzzy::add_pos_hd_cell ( Cell::PositiveHighDim * cll, const tag::MeshIsBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) == cll->meshes_same_dim.end() );
	
  make_deep_connections ( cll, this, tag::mesh_is_bdry );
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		assert ( face_p->cell_behind_within.find(this) ==
		         face_p->cell_behind_within.end()         );
		face_p->cell_behind_within[this] = cll;              }                          }


void Mesh::Fuzzy::remove_pos_hd_cell ( Cell::PositiveHighDim * cll, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) != cll->meshes_same_dim.end() );

	break_deep_connections ( cll, this, tag::mesh_is_not_bdry );
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		assert ( face_p->cell_behind_within.find(this) !=
		         face_p->cell_behind_within.end()         );
		assert ( face_p->cell_behind_within[this] == cll );
		face_p->cell_behind_within.erase(this);               }                          }

	
void Mesh::Fuzzy::remove_pos_hd_cell ( Cell::PositiveHighDim * cll, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) != cll->meshes_same_dim.end() );

	break_deep_connections ( cll, this, tag::mesh_is_bdry );
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		assert ( face_p->cell_behind_within.find(this) !=
		         face_p->cell_behind_within.end()         );
		assert ( face_p->cell_behind_within[this] == cll );
		face_p->cell_behind_within.erase(this);               }                          }

	
void Mesh::Fuzzy::add_neg_hd_cell ( Cell::NegativeHighDim * cll, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) == pos_cll->meshes_same_dim.end() );

	make_deep_connections_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry );
	
	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within.find(this) ==
		         rev_face->cell_behind_within.end() );
		rev_face->cell_behind_within[this] = cll;             }                           }	


void Mesh::Fuzzy::add_neg_hd_cell ( Cell::NegativeHighDim * cll, const tag::MeshIsBdry & )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) == pos_cll->meshes_same_dim.end() );

	make_deep_connections_rev ( cll, pos_cll, this, tag::mesh_is_bdry );
	
	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within.find(this) ==
		         rev_face->cell_behind_within.end() );
		rev_face->cell_behind_within[this] = cll;             }                           }	

	
void Mesh::Fuzzy::remove_neg_hd_cell ( Cell::NegativeHighDim * cll, const tag::MeshIsNotBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) != pos_cll->meshes_same_dim.end() );

	break_deep_connections_rev ( cll, pos_cll, this, tag::mesh_is_not_bdry );

	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within[this] == cll );
		rev_face->cell_behind_within.erase(this);              }                         }

	
void Mesh::Fuzzy::remove_neg_hd_cell ( Cell::NegativeHighDim * cll, const tag::MeshIsBdry & )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) != pos_cll->meshes_same_dim.end() );

	break_deep_connections_rev ( cll, pos_cll, this, tag::mesh_is_bdry );

	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell::Core * face_p = (*it).core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within[this] == cll );
		rev_face->cell_behind_within.erase(this);              }                         }

}  // anonymous namespace

	
//-----------------------------------------------------------------------------//




