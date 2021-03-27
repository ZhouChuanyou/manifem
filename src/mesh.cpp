
// mesh.cpp 2021.03.27

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

	
Cell::Core * const Cell::ghost { new Cell::Negative::Vertex ( tag::ghost ) };
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

bool Mesh::Sign::Positive::is_positive ( ) // static
{	return true;  }

bool Mesh::Sign::Negative::is_positive ( ) // NegativeMesh is a namespace
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

Cell::Core * Cell::Positive::Vertex::build_reverse ( )
// virtual from Cell::Core
{	this->reverse_p = new Cell::Negative::Vertex ( tag::reverse_of, this );
	return this->reverse_p;                                                    }

Cell::Core * Cell::Positive::Segment::build_reverse ( )
// virtual from Cell::Core
{	this->reverse_p = new Cell::Negative::Segment ( tag::reverse_of, this );
	return this->reverse_p;                                                     }

Cell::Core * Cell::Positive::HighDim::build_reverse ( )
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


size_t Cell::Positive::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Negative::Vertex::get_dim ( ) const  // virtual from Cell::Core
{	return 0;  }

size_t Cell::Positive::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Negative::Segment::get_dim ( ) const  // virtual from Cell::Core
{	return 1;  }

size_t Cell::Positive::HighDim::get_dim ( ) const  // virtual from Cell::Core
{	assert ( this->boundary_p );
	return this->boundary_p->get_dim_plus_one();  }

size_t Cell::Negative::HighDim::get_dim ( ) const  // virtual from Cell::Core
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
	
size_t Mesh::ZeroDim::number_of ( const tag::CellsOfDim &, size_t d )
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
	
size_t Mesh::Connected::OneDim::number_of ( const tag::CellsOfDim &, size_t d )
// virtual from Mesh::Core
{	if ( d == 1 ) return this->nb_of_segments;
	assert ( d == 0 );
	if ( first_ver == last_ver )  // closed loop
		return this->nb_of_segments;
	return this->nb_of_segments + 1;             }
	
size_t Mesh::Connected::HighDim::number_of ( const tag::CellsOfDim &, size_t d )
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
	
size_t Mesh::Fuzzy::number_of ( const tag::CellsOfDim &, size_t d )
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


Mesh Cell::Positive::Vertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                                                 }

Mesh Cell::Negative::Vertex::boundary ( )  // virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "vertices have no boundary" << std::endl;
	exit ( 1 );                                                                 }

Mesh Cell::Positive::Segment::boundary ( )  // virtual from Cell::Core
{	return Mesh ( tag::boundary_of, tag::positive, tag::segment, seg );  }

Mesh Cell::Negative::Segment::boundary ( )  // virtual from Cell::Core
{	assert ( this->reverse_p );
	Cell::Positive::Segment * rev_seg = Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( this->reverse_p );
	return Mesh ( tag::boundary_of, tag::positive, tag::segment, rev_seg,
                tag::reversed, tag::cells_surely_exist                  );  }

Mesh Cell::Positive::HighDim::boundary ( )  // virtual from Cell::Core
{	return Mesh ( tag::whose_core_is, this->boundary_p, tag::is_positive );   }

Mesh Cell::Negative::HighDim::boundary ( )  // virtual from Cell::Core
{	Cell::Positive::HighDim * cll = Util::assert_cast
		< Cell::Core*, Cell::Positive::HighDim* > ( this->reverse_p );
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


bool Cell::Positive::HighDim::belongs_to ( Mesh::Core * msh, const tag::Oriented & ) const
// virtual from Cell::Core

{	assert ( msh->get_dim_plus_one() == this->get_dim() + 1 );
	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap = this->meshes[0];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	if ( it == mmap.end() ) return false;
	const Cell::field_to_meshes & field = it->second;
	if ( field.counter_pos != 1 ) return false;
	assert ( field.counter_neg == 0 );
	return true;                                                                    }


bool Cell::Positive::HighDim::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	const std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		this->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::const_iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                     }


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
	return true;                                                                               }


bool Cell::Negative::belongs_to ( Mesh::Core * msh, const tag::NotOriented & ) const
// virtual from Cell::Core

{	assert ( this->reverse_p );
	assert ( this->reverse_p->is_positive() );
	Cell::Positive * rev = static_cast < Cell::Positive * > ( this->reverse_p );
	std::map < Mesh::Core *, Cell::field_to_meshes > & mmap =
		rev->meshes [ Mesh::diff ( msh->get_dim_plus_one(), this->get_dim() + 1 ) ];
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator it = mmap.find(msh);
	return ( it != mmap.end() );                                                                }

//-----------------------------------------------------------------------------//


Cell::Core * Cell::Core::tip ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::Core::base ()  // virtual
{	std::cout << __FILE__ << ":" <<__LINE__ << ": " << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Only segments have extremities." << std::endl;
	exit ( 1 );                                                                                   }

Cell::Core * Cell::Positive::Segment::tip () { return this->tip_p;  }
// virtual from Cell::Core, overridden here

Cell::Core * Cell::Positive::Segment::base () { return this->base_p;  }
// virtual from Cell::Core, overridden here

Cell::Core * Cell::Negative::Segment::tip ()  // virtual from Cell::Core, overridden here
{	assert ( this->reverse_p );
	Cell::Positive::Segment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( this->reverse_p );
	assert ( pos_seg->base_p );
	// assert ( pos_seg->base_p->reverse_p );
	return pos_seg->base_p->reverse_p;                                               }

Cell::Core * Cell::Negative::Segment::base ()  // virtual from Cell::Core, overridden here
{	assert ( this->reverse_p );
	Cell::Positive::Segment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( this->reverse_p );
	assert ( pos_seg->tip_p );
	// assert ( pos_seg->tip_p->reverse_p );
	return pos_seg->tip_p->reverse_p;                                                }

//-----------------------------------------------------------------------------//

#ifndef NDEBUG


std::string Cell::Positive::HighDim::get_name ()  // virtual from Cell::Core
{	return this->name;  }

std::string Cell::Negative::get_name ()  // virtual from Cell::Core
{	return "r" + this->reverse_p->name;  }


void Cell::Positive::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Positive::Vertex " << this->name << std::endl;
	if ( this->meshes.size() > 0 )
	{	if ( this->meshes[0].size() > 0 )
				std::cout << "meshes of index 0, dim 0 (segments disguised as meshes)" << std::endl;
		std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
			it = this->meshes[0].begin(), it_e = this->meshes[0].end();
		for ( ; it != it_e; it++ )
		{	Cell::Positive::Segment * seg = (Cell::Positive::Segment*) it->first;
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


void Cell::Negative::Vertex::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Cell::Negative::Vertex, reverse of "
            << this->reverse_p->name << std::endl;            }


void Cell::Positive::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Positive::Segment " << this->name << std::endl;
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


void Cell::Negative::Segment::print_everything ( )  // virtual from Cell::Core

{	std::cout << "this is Negative::Segment " << this->get_name() << std::endl;     }


void Cell::Positive::HighDim::print_everything ( )  // virtual from Cell::Core
	
{	size_t dim = this->get_dim();
	std::cout << "this is Cell::Positive::HighDim of dim "
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


void Cell::Negative::HighDim::print_everything ( )  // virtual from Cell::Core
	
{	std::cout << "this is Cell::Negative::HighDim "
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


void Cell::Positive::Vertex::glue_on_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Segment::glue_on_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->add_to ( (Mesh::Core*) this );  // now there are zero-dim meshes !!
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	// 'add_to' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::add_to or Cell::Negative::Vertex::add_to)
	this->glue_common ( ver );            }
	

void Cell::Positive::HighDim::glue_on_my_bdry ( Cell::Core * face )
// virtual from Cell::Core

{	assert ( this->get_dim() == face->get_dim() + 1 );
	face->add_to ( this->boundary_p );
	// 'add_to' is virtual, so the computer will choose the right version
	this->glue_common ( face );                  }


void Cell::Negative::glue_on_my_bdry ( Cell::Core * cll )
// virtual from Cell::Core

{	assert ( cll->reverse_p );
	assert ( this->reverse_p );
	this->reverse_p->glue_on_my_bdry ( cll->reverse_p );  }


void Cell::Positive::Vertex::cut_from_my_bdry ( Cell::Core * )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Vertices have no boundary." << std::endl;
	exit ( 1 );                                               }


void Cell::Positive::Segment::cut_from_my_bdry ( Cell::Core * ver )
// virtual from Cell::Core

{	assert ( ver->get_dim() == 0 );
	ver->remove_from ( (Mesh::Core*) this );  // now there are zero-dim meshes !!
	// ver.meshes[0] contains segments disguised as zero-dimensional meshes
	// 'remove_from' is virtual, so the computer will choose the right version
	// (Cell::Positive::Vertex::remove_from or Cell::Negative::Vertex::remove_from)
	this->cut_common ( ver );                           }
	

void Cell::Positive::HighDim::cut_from_my_bdry ( Cell::Core * face )
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


void Cell::Core::add_to_seg ( Cell::Positive::Segment * seg )
// virtual, overriden by Cell::***tive::Vertex
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a cell to a segment, use add_to_mesh." << std::endl;
	exit ( 1 );                                                                   }


void Cell::Core::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual, overriden by Cell::***tive::Vertex
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a cell from a segment, use remove_from_mesh." << std::endl;
	exit ( 1 );                                                                             }


void Cell::Positive::Vertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core, here overriden
	
{	assert ( seg );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	std::map < Cell::Positive::Segment *, short int > & tm0 = this->meshes_same_dim;
	assert ( tm0.find(seg) == tm0.end() );
	assert ( seg->tip_p == nullptr );
	seg->tip_p = this;
	tm0.emplace ( std::piecewise_construct,
		std::forward_as_tuple(seg), std::forward_as_tuple(1) );
	seg->deep_connections ( this, Mesh::action_add );                                  }


void Cell::Positive::Vertex::add_to_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                  }


void Cell::Positive::Vertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core, here overriden
	
{	assert ( seg );
	// assert that 'this' vertex belongs to the boundary of 'seg'
	std::map < Cell::Positive::Segment *, short int > & tm0 = this->meshes_same_dim;
	std::map<Cell::Positive::Segment*,short int>::iterator it = tm0.find ( seg );
	assert ( it != tm0.end() );
	assert ( seg->tip_p == this );
	seg->tip_p = nullptr;
	assert ( it->second == 1 );
	tm0.erase ( it );
	seg->deep_connections ( this, Mesh::action_remove );                             }


void Cell::Positive::Vertex::remove_from_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
            << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Negative::Vertex::add_to_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core, here overriden

{ assert ( seg );
	assert ( this->reverse_p );
	Cell::Positive::Vertex * pos_ver = Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_p );
	// assert that 'this' vertex does not belong yet to the boundary of 'seg'
	std::map < Cell::Positive::Segment *, short int > & pvm0 = pos_ver->meshes_same_dim;
	assert ( pvm0.find(seg) == pvm0.end() );
	assert ( seg->base_p == nullptr );
	seg->base_p = this;
	pvm0.emplace ( std::piecewise_construct,
		std::forward_as_tuple(seg), std::forward_as_tuple(-1) );
	seg->deep_connections ( pos_ver, Mesh::action_add_rev );                               }


void Cell::Negative::Vertex::add_to_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot add a vertex to a mesh, use add_to_seg." << std::endl;
	exit ( 1 );                                                                 }


void Cell::Negative::Vertex::remove_from_seg ( Cell::Positive::Segment * seg )
// virtual from Cell::Core, here overriden

{ assert ( seg );
	assert ( this->reverse_p );
	Cell::Positive::Vertex * pos_ver = Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( this->reverse_p );
	// assert that 'this' vertex belongs to the mesh 'msh'
	std::map < Cell::Positive::Segment *, short int > & pvm0 = pos_ver->meshes_same_dim;
	std::map<Cell::Positive::Segment*,short int>::iterator it = pvm0.find ( seg );
	assert ( it != pvm0.end() );
	assert ( seg->base_p );
	seg->base_p = nullptr;
	assert ( it->second == -1 );
	pvm0.erase ( it );
	seg->deep_connections ( pos_ver, Mesh::action_remove_rev );                      }


void Cell::Negative::Vertex::remove_from_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Cannot remove a vertex from a mesh, use remove_from_seg." << std::endl;
	exit ( 1 );                                                                            }


void Cell::Positive::Segment::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_seg ( this );  }


void Cell::Positive::Segment::remove_from_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_seg ( this );  }


void Cell::Negative::Segment::add_to_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_seg ( this );  }


void Cell::Negative::Segment::remove_from_mesh ( Mesh::Core * mmsh )
// virtual from Cell::Core

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_seg ( this );  }


void Cell::Positive::HighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Core

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_pos_hd_cell ( this );  }


void Cell::Positive::HighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Positive::NotVertex

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_pos_hd_cell ( this );  }


void Cell::Negative::HighDim::add_to_mesh ( Mesh::Core * msh )
// virtual from Cell::Negative::NotVertex

// add 'this' cell to the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->add_neg_hd_cell ( this );  }


void Cell::Negative::HighDim::remove_from_mesh ( Mesh::Core * msh )
// virtual from Cell::Negative::NotVertex

// remove 'this' cell from the mesh 'msh'
// if 'msh' is the boundary of a cell, use instead 'cut_from_bdry_of'

// there are many types of meshes, so we call a virtual method

{	assert ( msh ); msh->remove_neg_hd_cell ( this );  }

//-----------------------------------------------------------------------------//


void Mesh::ZeroDim::add_pos_seg ( Cell::Positive::Segment * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_pos_seg ( Cell::Positive::Segment * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_neg_seg ( Cell::Negative::Segment * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::remove_neg_seg ( Cell::Negative::Segment * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                             }
	
void Mesh::ZeroDim::add_pos_hd_cell ( Cell::Positive::HighDim * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::remove_pos_hd_cell ( Cell::Positive::HighDim * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::add_neg_hd_cell ( Cell::Negative::HighDim * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }
	
void Mesh::ZeroDim::remove_neg_hd_cell ( Cell::Negative::HighDim * )
// virtual from Mesh::Core, here execution forbidden
{	std::cout << __FILE__ << ":" <<__LINE__ << ": "
	          << __extension__ __PRETTY_FUNCTION__ << ": ";
	std::cout << "Something is very wrong, "
	          << "execution should never reach this point." << std::endl;
	exit ( 1 );                                                            }


void Mesh::Fuzzy::add_pos_seg ( Cell::Positive::Segment * seg )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' does not belong yet to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) == seg->meshes_same_dim.end() );

	Mesh::Fuzzy::deep_connections_add ( this, seg, seg );  // check !!

	assert ( seg->base_p->cell_behind_within.find(this) ==
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within[msh] = seg;
	assert ( seg->tip_p->cell_behind_within.find(this) ==
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within[this] = seg;                                 }


void Mesh::Fuzzy::remove_pos_seg ( Cell::Positive::Segment * seg )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->base_p );
	assert ( seg->tip_p );
	// assert that 'seg' belongs to 'this' mesh
	assert ( seg->meshes_same_dim.find(this) != seg->meshes_same_dim.end() );

	Mesh::Fuzzy::deep_connections_remove ( this, seg, seg );  // check !!

	assert ( seg->base_p->cell_behind_within.find(this) !=
	         seg->base_p->cell_behind_within.end()         );
	seg->base_p->cell_behind_within.erase(this);
	assert ( seg->tip_p->cell_behind_within.find(this) !=
	         seg->tip_p->cell_behind_within.end()         );
	seg->tip_p->cell_behind_within.erase(this);                                 }
	

void Mesh::Fuzzy::add_neg_seg ( Cell::Negative::Segment * seg )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::Positive::Segment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_p );
	// assert that 'pos_seg' does not belong yet to 'this' mesh
	assert ( seg_pos->meshes_same_dim.find(this) == seg_pos->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	Mesh::Fuzzy::deep_connections_add_rev ( this, pos_seg, seg );  // check !!

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	pos_seg->base_p->reverse_p->cell_behind_within[this] = seg;
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) ==
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	pos_seg->tip_p->reverse_p->cell_behind_within[this] = seg;                            }

	
void Mesh::Fuzzy::remove_neg_seg ( Cell::Negative::Segment * seg )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == 2 );
	assert ( seg->reverse_p );
	Cell::Positive::Segment * pos_seg = Util::assert_cast
		< Cell::Core*, Cell::Positive::Segment* > ( seg->reverse_p );
	// assert that 'pos_seg' belongs to 'this' mesh
	assert ( pos_seg->meshes_same_dim.find(this) != pos_seg->meshes_same_dim.end() );

	assert ( pos_seg->base_p );
	assert ( pos_seg->tip_p );
	assert ( pos_seg->base_p->reverse_p );
	assert ( pos_seg->tip_p->reverse_p );

	Mesh::Fuzzy::deep_connections_remove_rev ( this, pos_seg, seg );  // check !!

	assert ( pos_seg->base_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->base_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->base_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->base_p->reverse_p->cell_behind_within.erase(this);
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within.find(this) !=
	         pos_seg->tip_p->reverse_p->cell_behind_within.end()         );
	assert ( pos_seg->tip_p->reverse_p->cell_behind_within[this] == seg );
	pos_seg->tip_p->reverse_p->cell_behind_within.erase(this);                           }
	
	
void Mesh::Fuzzy::add_pos_hd_cell ( Cell::Positive::HighDim * cll )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' does not belong yet to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) == cll->meshes_same_dim.end() );
	
  Mesh::Fuzzy::deep_connections_add ( this, cll, cll );  // check !!
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell face = *it;
		Cell::Core * face_p = face.core;
		assert ( face_p->cell_behind_within.find(this) ==
		         face_p->cell_behind_within.end()         );
		face_p->cell_behind_within[this] = cll;              }                          }


void Mesh::Fuzzy::remove_pos_hd_cell ( Cell::Positive::HighDim * cll )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	// assert that 'cll' belongs to 'this' mesh
	assert ( cll->meshes_same_dim.find(this) != cll->meshes_same_dim.end() );

	Mesh::Fuzzy::deep_connections_remove ( this, cll, cll );  // check !!
	
	Mesh::Core * bdry = cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell face = *it;
		Cell::Core * face_p = face.core;
		assert ( face_p->cell_behind_within.find(this) !=
		         face_p->cell_behind_within.end()         );
		assert ( face_p->cell_behind_within[this] == cll );
		face_p->cell_behind_within.erase(this);               }                          }

	
void Mesh::Fuzzy::add_neg_hd_cell ( Cell::Negative::HighDim * cll )
// virtual from Mesh::Core

{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'pos_cll' does not belong yet to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) == pos_cll->meshes_same_dim.end() );

	Mesh::Fuzzy::deep_connections_add_rev ( this, pos_cll, cll );  // check !!
	
	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell face = *it;
		Cell::Core * face_p = face.core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within.find(this) ==
		         rev_face->cell_behind_within.end() );
		rev_face->cell_behind_within[this] = cll;             }                           }	

	
void Mesh::Fuzzy::remove_neg_hd_cell ( Cell::Negative::HighDim * cll )
// virtual from Mesh::Core
	
{	assert ( this->get_dim_plus_one() == cll->get_dim() + 1 );
	assert ( cll->reverse_p );
	Cell::Positive * pos_cll = Util::assert_cast
		< Cell::Core*, Cell::Positive* > ( cll->reverse_p );
	// assert that 'cll' belongs to 'this' mesh
	assert ( pos_cll->meshes_same_dim.find(this) != pos_cll->meshes_same_dim.end() );

	Mesh::Fuzzy::deep_connections_remove_rev ( this, pos_cll, cll );  // check !!

	Mesh::Core * bdry = pos_cll->boundary_p;  assert ( bdry );
	assert ( bdry->get_dim_plus_one() == this->get_dim() );
	CellIterator it = bdry->iterator ( tag::over_cells_of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell face = *it;
		Cell::Core * face_p = face.core;
		Cell::Core * rev_face = face_p->reverse_p;
		assert ( rev_face->cell_behind_within[this] == cll );
		rev_face->cell_behind_within.erase(this);              }                         }

	
//-----------------------------------------------------------------------------//

