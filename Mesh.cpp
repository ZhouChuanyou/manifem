// src/manifem/Mesh.cpp 2019.08.21

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "assert.h"
#include "Mesh.h"

using namespace ManiFEM;

// This file is organized in sections.
// Sections correspond to concepts, not necessarily to classes.


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//         HHH    HHHHHHH   H       H        HHHH                           //
//        H   H   H         H       H       H    H                    H     //
//       H        H         H       H       H                         H     //
//       H        HHHHH     H       H        H            HH   HH    HH     //
//       H        H         H       H         HHH        H  H  H H  H H     //
//       H        H         H       H            H       H  H  H H  H H     //
//        H   H   H         H       H       H    H        HHH  H H   HH     //
//         HHH    HHHHHHH   HHHHHH  HHHHHH   HHHH                           //
//                                                                          //
//                                                                          //
//           H     H   HHHHHH    HHHH    H    H   HHHHHHH    HHHH           //
//           HH   HH   H        H    H   H    H   H         H    H          //
//           H H H H   H        H        H    H   H         H               //
//           H  H  H   HHHH      HHH     HHHHHH   HHHHH      HHH            //
//           H     H   H            H    H    H   H             H           //
//           H     H   H             H   H    H   H              H          //
//           H     H   H        H    H   H    H   H         H    H          //
//           H     H   HHHHHH    HHHH    H    H   HHHHHHH    HHHH           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////
////////////////    constructors    ////////////////////
////////////////////////////////////////////////////////


Cell::Cell ( short int d, bool is_segment )

// second argument defaults to false
// accepts d == 1 if second argument is true, but builds an inconsistent segment !
// so do not use directly, use factory function Cell::segment instead
// (see paragraph 5.6 of the manual)
// bottom line : never use the second argument (only used once in Segment::Segment)
	
:	dim ( d ),
	hidden_boundary ( NULL ),
	hidden_reverse ( NULL ),
	real_heap ( Cell::real_heap_size[d] ),
	size_t_heap ( Cell::size_t_heap_size[d] ),
	short_int_heap ( Cell::short_int_heap_size[d] ),
	positive ( this )

{	assert ( is_segment or ( d != 1 ) );
	std::list < void(*)(Cell*,void*) > * ini = this->init_cell[d];
	std::list<void(*)(Cell*,void*)>::iterator it_ini;
	std::list < void* > * dat = this->data_for_init[d];
	std::list<void*>::iterator it_dat;
	for ( it_ini = ini->begin(), it_dat = dat->begin();
				it_ini != ini->end(); it_ini++, it_dat++ )
	{	assert ( it_dat != dat->end() );
		(*it_ini)(this,*it_dat);         }
	assert ( it_dat == dat->end() );                                }


Cell::Cell ( const tag::ReverseOf &, Cell & direct_cell )

// for building the reverse cell

:	dim ( direct_cell.dim ),
	real_heap ( Cell::real_heap_size_r[direct_cell.dim] ),
	size_t_heap ( Cell::size_t_heap_size_r[direct_cell.dim] ),
	short_int_heap ( Cell::short_int_heap_size_r[direct_cell.dim] ),
	positive ( &direct_cell ), hidden_reverse ( &direct_cell )

{	assert ( direct_cell.is_positive () );
	direct_cell.hidden_reverse = this;
	if ( direct_cell.hidden_boundary == NULL )
		this->hidden_boundary = NULL;
	else
	{	this->hidden_boundary = & ( direct_cell.hidden_boundary->reverse() );
		this->hidden_boundary->cell_enclosed = this;                           }
	std::list < void(*)(Cell*,void*) > * ini = this->init_cell_r[this->dim];
	std::list<void(*)(Cell*,void*)>::iterator it_ini;
	std::list < void* > * dat = this->data_for_init_r[this->dim];
	std::list<void*>::iterator it_dat;
	for ( it_ini = ini->begin(), it_dat = dat->begin();
				it_ini != ini->end(); it_ini++, it_dat++ )
	{	assert ( it_dat != dat->end() );
		(*it_ini) ( this, *it_dat );       }
	assert ( it_dat == dat->end() );                             }
// shouldn't we give special treatment to segments ?


Mesh::Mesh ( const tag::ReverseOf &, Mesh & direct_mesh )

// for building the reverse mesh

:	dim ( direct_mesh.dim ), hidden_reverse ( &direct_mesh ), positive ( &direct_mesh )

// there is no reason to ever call this constructor directly
// it is invoked from Mesh::reverse(), which in turn is
// usually called from Cell::Cell(tag::ReverseOf,Cell&) above

{	assert ( direct_mesh.is_positive() );
	direct_mesh.hidden_reverse = this;

	// We want to make sure that all (reverse) cells exist,
	// and to create the 'cell_behind' attribute.
	// We cannot use an iterator over 'this' negative mesh because
	// the iterator assumes the reverse cells exist (uses 'hidden_reverse').
	// Instead, we use an iterator over the direct_mesh and then
	// we call the 'reverse' method on each cell.
	CellIterator it = direct_mesh.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell & rev_cll = (*it).reverse();
		// now we are sure the reverse cell exists
		if (this->dim == 0) continue;
		CellIterator itt = rev_cll.boundary().iter_over
			( tag::cells, tag::of_max_dim, tag::oriented );
		// is this really necessary ?
		for ( itt.reset(); itt.in_range(); itt++ )
		{	// optimizar !
			Cell & face = *itt;
			assert ( face.cell_behind_within.find(this) == face.cell_behind_within.end() );
			face.cell_behind_within[this] = & rev_cll;                                       }  }
}


////////////////////////////////////////////////////////
/////////////    factory functions    ////////////////
////////////////////////////////////////////////////////


Cell & Cell::point ()
	
{	return * ( new Cell ( 0 ) );  }

Cell & Cell::segment ( Cell & a, Cell & b )

{	assert ( ! ( a.is_positive() ) );
	assert ( b.is_positive() );
	// we build 'seg' as a Segment, but we return it as a Cell
	Cell & seg =  * ( (Cell*) new Segment () );
	// we build the boundary, initially empty :
	Mesh::build_boundary ( seg );
	// we add faces :
	a.glue_on_bdry_of ( seg );
	b.glue_on_bdry_of ( seg );
	return seg;                                                  }

Cell & Cell::triangle ( Cell & a, Cell & b, Cell & c )

{	Cell & tri = * ( new Cell ( 2 ) );
	// we build the boundary, initially empty :
	Mesh::build_boundary ( tri );
	// we add faces :
	a.glue_on_bdry_of ( tri );
	b.glue_on_bdry_of ( tri );
	c.glue_on_bdry_of ( tri );
	return tri;                                    }

Cell & Cell::rectangle ( Cell & south, Cell & east, Cell & north, Cell & west )

{	Cell & rec = * ( new Cell ( 2 ) );
	// we build the boundary, initially empty :
	Mesh::build_boundary ( rec );
	// we add faces :
	south.glue_on_bdry_of ( rec);
	east.glue_on_bdry_of ( rec );
	north.glue_on_bdry_of ( rec );
	west.glue_on_bdry_of ( rec );
	return rec;                                    }


#ifndef NDEBUG

Cell & Cell::point (std::string s)

{	Cell & p = Cell::point ();
	p.hidden_name = s;
	return p;                   }

Cell & Cell::segment ( std::string s, Cell & a, Cell & b )

{	Cell & seg = Cell::segment ( a, b );
	seg.hidden_name = s;
	return seg;                        }

Cell & Cell::triangle ( std::string s, Cell & a, Cell & b, Cell & c )

{	Cell & tri = Cell::triangle ( a, b, c );
	tri.hidden_name = s;
	return tri;                             }

Cell & Cell::rectangle ( std::string s, Cell & south, Cell & east, Cell & north, Cell & west )

{	Cell & rec = Cell::rectangle ( south, east, north, west );
	rec.hidden_name = s;
	return rec;                                              }

#endif


////////////////////////////////////////////////////////////
/////////////////////  more stuff  /////////////////////////
////////////////////////////////////////////////////////////


// The user should set the correct value of Mesh::intended_dimension at the beginning
// of the main program (see comments about Mesh::initialize below).
// We initialize it with value 3, however, just in case the user forgets about it.
// Note : this is the topological dimension (may be less than the number of coordinates).

Cell * Mesh::current_cell_for_prescribe;

hidden::VectorOfIntegers < size_t > Cell::real_heap_size,	Cell::size_t_heap_size,
	Cell::short_int_heap_size, Cell::real_heap_size_r,	Cell::size_t_heap_size_r,
	Cell::short_int_heap_size_r;

hidden::IndexedList < void(*)(Cell*,void*) > Cell::init_cell, Cell::init_cell_r;
hidden::IndexedList < void* > Cell::data_for_init, Cell::data_for_init_r;

Cell * FiniteElement::current_cell;
Manifold * Mesh::environment;

FunctionOnMesh::Function
FunctionOnMesh::zero ( * ( (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Constant(0.) ) ),
FunctionOnMesh::one ( * ( (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Constant(1.) ) ),
FunctionOnMesh::minus_one ( * ( (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Constant(-1.) ) );
// std::map < std::string, FE_function* > FE_function::variables;

VariationalProblem * VariationalProblem::current_var_pb;
Mesh * VariationalProblem::current_mesh;


////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NDEBUG

// Cell::name, Mesh::name, Cell::print_meshes_above, Mesh::print_cells,
// Cell::print_everything and Mesh::print_everything only exist for debugging purposes

std::string Cell::name ()
{	if ( this->is_positive() ) return this->hidden_name;
	else return this->reverse().hidden_name + "-";       }

std::string Mesh::name ()
{	if (this->cell_enclosed == NULL) return "no-name";
	else return "bdry of " + this->cell_enclosed->name();  }


void Cell::print_meshes_above ()

{	assert ( this->is_positive() );
	std::cout << "as malhas acima de " << this->name() << " : " << std::endl;
	hidden::IndexedMap<Mesh*,Cell::field_to_meshes>::iterator dim_iter = this->meshes.begin();
	for ( ; dim_iter != this->meshes.end(); ++dim_iter)
	{	short int dim = dim_iter->first;
		std::cout << "dimensão " << dim << std::endl;
		std::map<Mesh*,Cell::field_to_meshes>::iterator map_iter = this->meshes[dim]->begin();
		for (; map_iter != this->meshes[dim]->end(); ++map_iter)
			std::cout  << map_iter->first->name ()
			  << " " << map_iter->second.counter_pos << " vezes; "
			  << map_iter->first->name () << " reversed "
			  << map_iter->second.counter_neg << " vezes" << std::endl << std::flush;                  }
} // end of Cell::print_meshes_above


void Mesh::print_cells ()

{	assert ( this->is_positive() );
	std::cout << "as células de " << this->name() << " : " << std::endl;
	hidden::IndexedList<Cell*>::reverse_iterator map_iter = this->cells.rbegin();
	for ( ; map_iter != this->cells.rend(); ++map_iter)
	{	short int dim = map_iter->first;
		std::cout << "dimensão " << dim << std::endl;
		// usar it = this->iter_over ("cells of dimension", dim, "not oriented");
		std::list<Cell*>::iterator list_iter = this->cells[dim]->begin();
		for ( ; list_iter != this->cells[dim]->end(); ++list_iter )
		{	std::cout << (*list_iter)->name() << " "
		        	<< (*(*list_iter)->meshes[this->dim])[this].counter_pos
				<< " vezes" << ", " << (*list_iter)->name() << "- "
				<< (*(*list_iter)->meshes[this->dim])[this].counter_neg
				<< " vezes, " << std::endl;
			if ( dim == 1 )
			{	std::cout << "base : " << (*list_iter)->base().name() << std::endl;
				std::cout << "tip : " <<  (*list_iter)->tip().name() << std::endl;   }   }  }
} // end of Mesh::print_cells


void Cell::print_everything ()

{	std::cout << "cell " << this->name() << std::endl;
	if ( this->hidden_reverse != NULL )
		std::cout << "reverse : " << this->hidden_reverse->name() << std::endl;
	if ( this->is_positive() )
	{	hidden::IndexedMap <Mesh*, Cell::field_to_meshes> :: iterator
			dim_iter = this->meshes.begin(), dim_e = this->meshes.end();
		for ( ; dim_iter != dim_e; dim_iter++ )
		{	std::cout << "meshes of dimension " << dim_iter->first << std::endl;
			std::map <Mesh*, Cell::field_to_meshes> :: iterator
				map_iter = dim_iter->second->begin(),
				map_e = dim_iter->second->end();
			for ( ; map_iter != map_e; map_iter++ )
			{	Mesh *msh = map_iter->first;
				Cell::field_to_meshes &field = map_iter->second;
				std::cout << msh->name() << " : " << field.counter_pos << ", "
					<< field.counter_neg << ", 'where' points towards "
					<< (*field.where)->name() << std::endl;                   }  }  }
	if (this->cell_behind_within.size() > 0)
	{	std::map <Mesh*, Cell*> :: iterator
			it = this->cell_behind_within.begin(),
			e = this->cell_behind_within.end();
		for ( ; it != e; it++ )
			std::cout << "cell_behind, within mesh "
				<< it->first->name() << " : " << it->second->name() << std::endl;  }
	if (this->dim == 1) // base and tip of segment
		std::cout << "base : " << this->base().name() << ", tip : " << this->tip().name() << std::endl;
} // end of Cell::print_everything


void Mesh::print_everything ()

{	std::cout << "mesh " << this->name() << std::endl;
	if ( this->hidden_reverse != NULL )
		std::cout << "reverse : " << this->hidden_reverse->name() << std::endl;
	if ( this->cell_enclosed != NULL )
		std::cout << "cell enclosed : " << this->cell_enclosed->name() << std::endl;
	if ( ! ( this->is_positive() ) ) return;
	hidden::IndexedList <Cell*> :: reverse_iterator dim_iter = this->cells.rbegin(),
		dim_e = this->cells.rend();
	for ( ; dim_iter != dim_e; dim_iter++ )
	{	std::cout << "cells of dimension " << dim_iter->first << " : ";
		std::list <Cell*> :: iterator list_iter = dim_iter->second->begin(),
			list_e = dim_iter->second->end();
		for ( ; list_iter != list_e; list_iter++ )
			std::cout << (*list_iter)->name() << " ";
		std::cout << std::endl;                                                           }
} // end of Mesh::print_everything

#endif


///////////////////////////////////////////////////////////////////
////////////    add/remove a cell to/from a mesh    ///////////////
///////////////////////////////////////////////////////////////////
////////////     this is where the magic happens     //////////////
//////////// (and also where we get our hands dirty) //////////////
///////////////////////////////////////////////////////////////////


void Cell::glue_on_bdry_of ( Cell & cll )

// glue 'this' face on the boundary of cell 'cll'

{	assert ( cll.hidden_boundary != NULL );
	assert ( cll.dim == this->dim + 1 );

	this->add_to ( cll.boundary() );

	// perhaps we should do this only to 'this->positive'
	this->glue_on_bdry_core ( cll );
	if ( cll.hidden_reverse != NULL )
		this->reverse().glue_on_bdry_core ( * ( cll.hidden_reverse ) );   }


void Cell::glue_on_bdry_core ( Cell & cll )

// used only in Cell::glue_on_bdry_of

{	short int cll_d = cll.dim;
	bool cll_pos = cll.is_positive();
	std::map <Mesh*, Cell::field_to_meshes> * cmd;
	if ( cll_pos ) cmd = cll.meshes[cll_d];
	else // when cll is negative
	{	assert ( cll.hidden_reverse != NULL );
		cmd = cll.hidden_reverse->meshes[cll_d];  }
	std::map <Mesh*, Cell::field_to_meshes> :: iterator it = cmd->begin(), e = cmd->end();
	for ( ; it != e; ++it)
	{	Mesh * msh_p = it->first;
		assert ( msh_p->is_positive() );
		Cell::field_to_meshes & field = it->second;
		bool reverse;
		if ( field.counter_pos != 0 )
		{	assert ( field.counter_pos == 1 );
			assert ( field.counter_neg == 0 );
			reverse = not cll_pos;              }
		else
		{	assert ( field.counter_pos == 0 );
			assert ( field.counter_neg == 1 );
			reverse = cll_pos;                  }
		if ( reverse )
		{	assert ( this->hidden_reverse != NULL );
			assert ( cll.hidden_reverse != NULL );
			assert ( this->hidden_reverse->cell_behind_within.find(msh_p) ==
		           this->hidden_reverse->cell_behind_within.end()         );
			this->hidden_reverse->cell_behind_within[msh_p] = cll.hidden_reverse;  }
		else
		{	assert ( this->cell_behind_within.find(msh_p) ==
		           this->cell_behind_within.end()         );
			this->cell_behind_within[msh_p] = & cll;              }
	} // end of for
	
	// special treatment is given to segments :
	if ( cll.dim == 1 )
	{	Segment & seg = * ( (Segment*) (& cll) );
		if ( this->is_positive() ) seg.hidden_tip = this;
		else seg.hidden_base = this;                        }
			
} // end of Cell::glue_on_bdry_core


void Cell::cut_from_bdry_of ( Cell & cll )

// cut 'this' face from the boundary of cell 'cll'

{	assert ( cll.hidden_boundary != NULL );
	assert ( cll.dim == this->dim + 1 );

	this->remove_from ( * ( cll.hidden_boundary ) );
	
	this->cut_from_bdry_core ( cll );
	if ( cll.hidden_reverse != NULL )
		this->reverse().cut_from_bdry_core ( * ( cll.hidden_reverse ) );    }


void Cell::cut_from_bdry_core ( Cell & cll )

// used only in Cell::cut_from_bdry_of

{	short int cll_d = cll.dim;
	bool cll_pos = cll.is_positive();
	std::map <Mesh*, Cell::field_to_meshes> *cmd;
	if ( cll_pos ) cmd = cll.meshes[cll_d];
	else // when cll is negative
	{	assert ( cll.hidden_reverse != NULL );
		cmd = cll.hidden_reverse->meshes[cll_d];  }
	std::map <Mesh*, Cell::field_to_meshes> :: iterator it = cmd->begin(), e = cmd->end();
	for ( ; it != e; ++it)
	{	Mesh *msh = it->first;
		assert ( msh->is_positive() );
		Cell::field_to_meshes &field = it->second;
		bool reverse;
		if (field.counter_pos != 0)
		{	assert (field.counter_pos == 1);
			assert (field.counter_neg == 0);
			reverse = !cll_pos;                 }
		else
		{	assert (field.counter_pos == 0);
			assert (field.counter_neg == 1);
			reverse = cll_pos;                 }
		if (reverse)
		{	assert ( this->hidden_reverse != NULL );
			assert ( cll.hidden_reverse != NULL );
			assert ( this->hidden_reverse->cell_behind_within.find(msh) !=
		           this->hidden_reverse->cell_behind_within.end()         );
			this->hidden_reverse->cell_behind_within.erase (msh);               }
		else
		{	assert ( this->cell_behind_within.find(msh) !=
		           this->cell_behind_within.end()         );
			this->cell_behind_within.erase (msh);               }
	} // end of for
	
	// special treatment is given to segments :
	if ( cll.dim == 1 )
	{	Segment & seg = * ( (Segment*) (&cll) );
		if ( this->is_positive() ) seg.hidden_tip = NULL;
		else seg.hidden_base = NULL;                       }
} // end of Cell::cut_from_bdry_core


void Cell::add_to ( Mesh & msh )

// add 'this' cell to the mesh 'msh'
// four possible orientation combinations are taken into account
// (two directly and two through recursivity)

{	assert ( this->dim == msh.dim );
	if ( ! ( msh.is_positive() ) )
	{	assert ( msh.hidden_reverse != NULL );
		this->reverse().add_to ( * ( msh.hidden_reverse ) );
		return;                                         }

	if ( this->is_positive() )
	{	// assert that 'this' cell does not belong yet to the mesh 'msh'
		assert ( this->meshes[msh.dim]->find(&msh) == this->meshes[msh.dim]->end() );
		// optimizar !
		msh.deep_connections ( *this, Mesh::action_add );                               }
	else // when 'this' is negative
	{	// assert that the reverse of 'this' cell does not belong yet to the mesh 'msh'
		assert ( this->hidden_reverse->meshes[msh.dim]->find(&msh)
			== this->hidden_reverse->meshes[msh.dim]->end() );
		// optimizar !
		msh.deep_connections ( * ( this->hidden_reverse ), Mesh::action_add_rev );    }

	if ( this->dim > 0 )
	{	Mesh & bdry = this->boundary();
		// 'bdry' may be a negative mesh, so we use a MeshIterator.
		CellIterator it = bdry.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( it.reset(); it.in_range(); it.advance() )
		{	Cell & face = *it;
			assert ( face.cell_behind_within.find(&msh) == face.cell_behind_within.end() );
			face.cell_behind_within[&msh] = this;  	                                         }   }

} // end of Cell::add_to


void Cell::remove_from (Mesh & msh)

// remove 'this' cell from the mesh 'msh'
// four possible orientation combinations are taken into account
// (two directly and two through recursivity)

{	assert ( this->dim == msh.dim );
	if ( ! ( msh.is_positive() ) )
	{	assert ( msh.hidden_reverse != NULL );
		this->reverse().remove_from ( * ( msh.hidden_reverse ) );
		return;                                        }

	if ( this->is_positive() )
	{	// assert that 'this' cell belongs to the mesh 'msh' exactly once
#ifndef NDEBUG
		std::map <Mesh*, Cell::field_to_meshes>::iterator iter_find = 
			this->meshes[msh.dim]->find(&msh);
		assert ( iter_find != this->meshes[msh.dim]->end() ); 
		assert ( iter_find->second.counter_pos == 1 );
		assert ( iter_find->second.counter_neg == 0 );
#endif
		msh.deep_connections ( *this, Mesh::action_remove );  }
	else
	{	// assert that the reverse of 'this' belongs to the mesh 'msh' exactly once
#ifndef NDEBUG
		std::map <Mesh*, Cell::field_to_meshes>::iterator iter_find = 
			this->hidden_reverse->meshes[msh.dim]->find(&msh);
		assert ( iter_find != this->hidden_reverse->meshes[msh.dim]->end() ); 
		assert ( iter_find->second.counter_pos == 0 );
		assert ( iter_find->second.counter_neg == 1 );
#endif
		msh.deep_connections ( * ( this->hidden_reverse ), Mesh::action_remove_rev );  }
	
	if ( this->dim != 0 )
	{	Mesh & bdry = this->boundary();
		// 'bdry' may be a negative mesh, so we use a CellIterator.
		CellIterator it = bdry.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( it.reset(); it.in_range(); it.advance() )
		{	Cell & face = *it;
			assert ( face.cell_behind_within.find(&msh) != face.cell_behind_within.end() );
			face.cell_behind_within.erase ( &msh );                                          }   }

} // end of Cell::remove_from


void Mesh::action_add ( Cell & cll, Mesh & msh, short int cp, short int cn ) // static

// This function "makes the link" between a cell and a mesh.
// This "link" asserts that the cell belongs to the mesh.
// If the dimensions are equal, then there will be only one instance of the cell
// in the mesh, so either cp==1 and cn==0 or cp==0 and cn==1, depending on the
// orientations.
// If 'cll' has lower dimension, it may appear several times in the mesh 'msh',
// with both orientations. The counter 'cp' asserts how many times 'cll' appears
// in the mesh with direct orientation (the same as the mesh), while 'cn' says
// the same thing for the reverse orientation (opposite to the mesh).
// This function is used only in add_cell, passed to deep_connections.

{	assert ( cll.is_positive() ); assert ( msh.is_positive() );
//////////////////////////////////////////////////////////////////////////
	// inspired in item 24 of the book : Scott Meyers, Effective STL      //
	typedef std::map <Mesh*, Cell::field_to_meshes> maptype;              //
	maptype* cmd = cll.meshes[msh.dim];                                   //
	maptype::iterator lb = cmd->lower_bound(&msh);                        //
	if ( ( lb == cmd->end() ) || ( cmd->key_comp()(&msh,lb->first) ) )    //
	{	std::list <Cell*>* mcd = msh.cells[cll.dim];                        //
		mcd->push_front ( & cll );                                          //
		Cell::field_to_meshes field;                                        //
		field.counter_pos = cp;                                             //
		field.counter_neg = cn;                                             //
		field.where = mcd->begin();                                         //
		maptype::value_type v ( & msh, field );                             //
		cmd->insert(lb,v);                                                } //
	else                                                                  //
	{	lb->second.counter_pos += cp;                                       //
		lb->second.counter_neg += cn;   }                                   //
////////// code below is functionally equivalent to the above ////////////
//	if (cll.meshes[msh.dim]->find(&msh)==cll.meshes[msh.dim]->end())    //
//	{	msh.cells[cll.dim]->push_front(&cll);                             //
//		Cell::field_to_meshes field;                                      //
//		field.counter_pos = cp;                                           //
//		field.counter_neg = cn;                                           //
//		field.where = msh.cells[cll.dim]->begin();                        //
//		(*cll.meshes[msh.dim])[&msh] = field;        }                    //
//	else                                                                //
//	{	(*cll.meshes[msh.dim])[&msh].counter_pos += cp;                   //
//		(*cll.meshes[msh.dim])[&msh].counter_neg += cn;    }              //
//////////////////////////////////////////////////////////////////////////

} // end of Mesh::action_add


void Mesh::action_remove ( Cell & cll, Mesh & msh, short int cp, short int cn ) // static
// used only in remove_cell, passed to deep_connections

{	typedef std::map <Mesh*, Cell::field_to_meshes> maptype;
	maptype* cmd = cll.meshes[msh.dim];
	maptype::iterator cmdm = cmd->find(&msh);
	assert (cmdm != cmd->end());
	short int c_p = cmdm->second.counter_pos -= cp;
	short int c_n = cmdm->second.counter_neg -= cn;
	assert ( (c_p >= 0) && (c_n >= 0) );
	if ( (c_p == 0) && (c_n == 0) )
	{	std::list <Cell*>::iterator w = cmdm->second.where;
		msh.cells[cll.dim]->erase(w);
		cmd->erase(&msh);                                }
} // end of Mesh::action_remove


void Mesh::deep_connections ( Cell & cell,
	void (*action)(Cell&, Mesh&, short int, short int) )

// make or destroy connections when addin or removing a cell,
// according to the 'action' argument

{	// We build two lists, a list of cells "below" 'cell'
	// (that is, belonging to the boundary of 'cell')
	// and a list of meshes "above" 'this' mesh
	// (that is, meshes to which 'this->cell_enclosed' belongs).
	// We then use these two lists to create or destroy all connections.
	
	struct triplet_cell
	{	Cell *obj;
		short int counter_pos;
		short int counter_neg;     };
	std::list < triplet_cell* > all_cells;
	{ // just a block of code, for hiding variables
	triplet_cell *t = new triplet_cell;
	t->obj = & cell;
	t->counter_pos = 1;
	t->counter_neg = 0;
	all_cells.push_back(t);
	Mesh *cell_bdry = cell.hidden_boundary;
	if ( cell_bdry != NULL )
	{	hidden::IndexedList<Cell*>::iterator dim_iter = cell_bdry->cells.begin(),
			cbc_end = cell_bdry->cells.end();
		// we loop over all cells of cell_bdry (of all dimensions)	
		for ( ; dim_iter != cbc_end; ++dim_iter)
		{	short int d = dim_iter->first;
			std::list<Cell*>* cbcd = cell_bdry->cells[d];
			std::list<Cell*>::iterator list_iter = cbcd->begin(),
				cbcd_end = cbcd->end();
			// we now loop over all cells of given dimension 'd'
			for (; list_iter != cbcd_end; ++list_iter)
			{	Cell* lower_cell = *list_iter;
				std::map <Mesh*, Cell::field_to_meshes>*
					f = lower_cell->meshes[cell_bdry->dim];
				Cell::field_to_meshes fcb = (*f)[cell_bdry];
				triplet_cell *t = new triplet_cell;
				t->obj = lower_cell;
				t->counter_pos = fcb.counter_pos;
				t->counter_neg = fcb.counter_neg;
				all_cells.push_back(t);                         }     }             }
	} // end of block of code which hides variables
	
	struct triplet_mesh
	{	Mesh *obj;
		short int counter_pos;
		short int counter_neg;     };
	std::list < triplet_mesh* > all_meshes;
	{ // just a block of code, for hiding variables
	triplet_mesh *t = new triplet_mesh;
	t->obj = this;
	t->counter_pos = 1;
	t->counter_neg = 0;
	all_meshes.push_back(t);
	if ( this->cell_enclosed != NULL )
	{	hidden::IndexedMap<Mesh*,Cell::field_to_meshes>::iterator 
			 dim_iter = this->cell_enclosed->meshes.begin(),
			 cem_end = this->cell_enclosed->meshes.end();
		// we now loop over all meshes above 'this->cell_enclosed'
		for ( ; dim_iter != cem_end; ++dim_iter)
		{	short int d = dim_iter->first;
			std::map<Mesh*, Cell::field_to_meshes>*
				cemd = this->cell_enclosed->meshes[d];
			std::map<Mesh*,Cell::field_to_meshes>::iterator
				map_iter = cemd->begin(), cemd_end = cemd->end();
			// we now loop over all meshes of given dimension 'd'
			for (; map_iter != cemd_end; ++map_iter)
			{	triplet_mesh *t = new triplet_mesh;
				t->obj = map_iter->first;
				t->counter_pos = map_iter->second.counter_pos;
				t->counter_neg = map_iter->second.counter_neg;
				all_meshes.push_back(t);                         }   }    }
	} // end of block of code which hides variables

	std::list < triplet_cell* >::iterator
		cells_iter = all_cells.begin(),
		cells_end = all_cells.end();
	// we loop over cells "below" 'cell'
	for ( ; cells_iter != cells_end; ++cells_iter)
	{	std::list < triplet_mesh* >::iterator
			meshes_iter = all_meshes.begin(),
			meshes_end = all_meshes.end();
		// we loop over "superior" dimensions	
		for (; meshes_iter != meshes_end; ++meshes_iter)
		{	triplet_cell *lower_cell_tr = *cells_iter;
			triplet_mesh *higher_mesh_tr = *meshes_iter;
			short int ccp = lower_cell_tr->counter_pos;
			short int ccn = lower_cell_tr->counter_neg;
			short int mcp = higher_mesh_tr->counter_pos;
			short int mcn = higher_mesh_tr->counter_neg;
			action ( *(lower_cell_tr->obj), *(higher_mesh_tr->obj),
			         ccp*mcp + ccn*mcn, ccp*mcn + ccn*mcp           );  }  }

	// so many new triplets ! should we clean up after us ?
	// or declare them local within all_meshes and all_cells ?
	// thus, in the end, when these lists are destroyed, so are all triplets

} // end of Mesh::deep_connections


////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//     H  HHHHHHHHH  HHHHHHH   HHHH        H    HHHHHHHHH   HHH     HHHH       HHHH       //
//     H      H      H         H   HH     HHH       H      H   H    H   HH    H    H      //
//     H      H      H         H    H    H   H      H     H     H   H     H   H           //
//     H      H      HHHH      HHHHH     H   H      H     H     H   HHHHHH     H          //
//     H      H      H         H H      HHHHHHH     H     H     H   H H         HH        //
//     H      H      H         H  H     H     H     H     H     H   H  H          H       //
//     H      H      H         H   H   H       H    H      H   H    H   H     H    H      //
//     H      H      HHHHHHH   H    H  H       H    H       HHH     H    H     HHHH       //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////


// Iterators are classified into two types.
// MeshIterators iterate over lists of meshes, e.g. over all meshes above a given cell.
// CellIterators iterate over lists of cells,
// e.g. over all cells of given dimension belonging to a given mesh,
// or over all cells of given dimension "above" another, given, cell.

// eliminar tudo com min_dim ?


void Cell::reset_mesh_iter_min_dim ( MeshIterator *it )  // static

{	Cell::data_for_iter_min_dim *data = (Cell::data_for_iter_min_dim*) it->data;
	Cell *cll = (Cell*) it->base;
	data->mesh_map = cll->meshes [cll->dim];
	data->iter = data->mesh_map->begin();
	it->skip_null ();                                                              }


void Cell::advance_mesh_iter_min_dim ( MeshIterator *it )  // static

{	Cell::data_for_iter_min_dim *data = (Cell::data_for_iter_min_dim*) it->data;
	data->iter++;
	it->skip_null ();                                                              }


bool Cell::valid_mesh_iter_min_dim ( MeshIterator *it )  // static

{	Cell::data_for_iter_min_dim *data = (Cell::data_for_iter_min_dim*) it->data;
	return ( data->iter != data->mesh_map->end() );                                 }


Mesh * Cell::deref_iter_pos_min_dim ( MeshIterator* it )  // static

{	Cell::data_for_iter_min_dim *data = (Cell::data_for_iter_min_dim*) it->data;
	std::pair <Mesh*, Cell::field_to_meshes> p = *data->iter;
	Mesh *msh = p.first;
	assert ( msh->is_positive() );
	Cell::field_to_meshes &field = p.second;
	if ( field.counter_pos != 0 )
	{	assert ( field.counter_pos == 1 );
		assert ( field.counter_neg == 0 );
		return msh;                          }
	else
	{	assert ( field.counter_pos == 0 );
		assert ( field.counter_neg == 1 );
		return msh->hidden_reverse;         }                                       }


Mesh * Cell::deref_iter_neg_min_dim ( MeshIterator* it )  // static

{	Cell::data_for_iter_min_dim *data = (Cell::data_for_iter_min_dim*) it->data;
	std::pair <Mesh*, Cell::field_to_meshes> p = *data->iter;
	Mesh *msh = p.first;
	assert ( msh->is_positive() );
	Cell::field_to_meshes &field = p.second;
	if ( field.counter_pos != 0 )
	{	assert ( field.counter_pos == 1 );
		assert ( field.counter_neg == 0 );
		return msh->hidden_reverse;            }
	else
	{	assert ( field.counter_pos == 0 );
		assert ( field.counter_neg == 1 );
		return msh;                          }                                         }


void Cell::reset_cell_iter_min_dim ( CellIterator *it, Cell * p )  // static

{	assert ( p == NULL );
	Cell::reset_mesh_iter_min_dim ( (MeshIterator*) it->data );  }


void Cell::advance_cell_iter_min_dim ( CellIterator *it )  // static

{	Cell::advance_mesh_iter_min_dim ( (MeshIterator*) it->data );  };


bool Cell::valid_cell_iter_min_dim ( CellIterator *it )  // static

{	return Cell::valid_mesh_iter_min_dim ( (MeshIterator*) it->data );  };


Cell * Cell::deref_cell_iter_min_dim ( CellIterator *it )  // static

{	MeshIterator *iter = (MeshIterator*) it->data;
	Cell *cll = (**iter).cell_enclosed;
	assert ( cll != NULL );
	return cll;                                     };


void Mesh::reset_iter_max_dim ( CellIterator *it, Cell *p )  // static

{	assert ( p == NULL );
	Mesh::data_for_iter_max_dim *data = (Mesh::data_for_iter_max_dim*) it->data;
	Mesh *msh = (Mesh*) it->base;
	data->cell_list = msh->cells [msh->dim];
	data->iter = data->cell_list->begin();                                         }


void Mesh::advance_iter_max_dim ( CellIterator *it )  // static

{	Mesh::data_for_iter_max_dim *data = (Mesh::data_for_iter_max_dim*) it->data;
	data->iter++;                                                                  }


bool Mesh::valid_iter_max_dim ( CellIterator *it )  // static

{	Mesh::data_for_iter_max_dim *data = (Mesh::data_for_iter_max_dim*) it->data;
	return ( data->iter != data->cell_list->end() );                               }


Cell * Mesh::deref_iter_pos_max_dim ( CellIterator *it )  // static
	
{	Mesh::data_for_iter_max_dim *data = (Mesh::data_for_iter_max_dim*) it->data;
	Cell *cll = *data->iter;
	assert ( cll->is_positive() );
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->is_positive() );
	Cell::field_to_meshes &field = (*cll->meshes[msh->dim])[msh];
	if ( field.counter_pos != 0 )
	{	assert ( field.counter_pos == 1 );
		assert ( field.counter_neg == 0 );
		return cll;                         }
	else
	{	assert ( field.counter_pos == 0 );
		assert ( field.counter_neg == 1 );
		assert ( cll->hidden_reverse != NULL );
		return cll->hidden_reverse;              }                                   }


Cell * Mesh::deref_iter_neg_max_dim ( CellIterator *it )  // static

{	Mesh::data_for_iter_max_dim *data = (Mesh::data_for_iter_max_dim*) it->data;
	Cell *cll = *data->iter;
	assert ( cll->is_positive() );
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->is_positive() );
	// return msh->orient_cell(cll)->reverse();
	Cell::field_to_meshes &field = (*cll->meshes[msh->dim])[msh];
	if ( field.counter_pos != 0 )
	{	assert ( field.counter_pos == 1 );
		assert ( field.counter_neg == 0 );
		assert ( cll->hidden_reverse != NULL );
		return cll->hidden_reverse;              }
	else
	{	assert ( field.counter_pos == 0 );
		assert ( field.counter_neg == 1 );
		return cll;                         }                                          }

////////////////////////////////////////////////////////////////////////////////////


void Mesh::reset_iter_along_ver ( CellIterator *it, Cell * p )  // static

// fazer ao contrario : chamar reset_ver do reset_seg
// e aqui, no reset_ver, procurar first point
	
{	Mesh::reset_iter_along_seg ( it, p );
	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
 	data->seg_got_out = ( data->current_seg == NULL );                             }


void Mesh::reset_iter_along_ver_rev ( CellIterator *it, Cell * p )  // static

{	Mesh::reset_iter_along_seg_rev ( it, p );
	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	data->seg_got_out = ( data->current_seg == NULL );                             }


void Mesh::reset_iter_along_seg ( CellIterator *it, Cell *p )  // static

// if p is null, initialize the iterator at the beginning of the curve
	
{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->cells[msh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	Mesh *omsh = data->oriented_base;
	std::string key;
	if ( p == NULL )
	// we try to use attribute "first point" or "last point"
	{	if ( omsh->is_positive() ) key = "first point";
		else key = "last point";
		std::map<std::string,void*>::iterator it = msh->hook.find(key);
		if ( it != msh->hook.end() )     // found
			p = (Cell*) it->second;                                        }
	if ( p == NULL )
	// p is still NULL, which means that the mesh has no "first point" attribute
	// (or "last point") thus we need to search for it
	// if the chain is not open (if it's a loop) this search may last forever !
	{	p = *(msh->cells[0]->begin());
	 	Cell *s = omsh->cell_behind ( p, tag::may_not_exist );
	 	while ( s )
	 	{	p = s->base().hidden_reverse;
	 		assert ( p != NULL );
	 		s = omsh->cell_behind( p, tag::may_not_exist );  }
		msh->hook[key] = (void*) p;                             }
	Mesh::reset_iter_around ( it, p );                                             }


void Mesh::reset_iter_along_seg_rev ( CellIterator *it, Cell *p )  // static

// if p is null, initialize the iterator at the beginning of the curve
{	Mesh::data_for_iter_contour *data =
		(Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->cells[msh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	Mesh *omsh = data->oriented_base;
	std::string key;
	if ( p == NULL )
	// we try to use attribute "first point" or "last point"
	{	if ( omsh->is_positive() ) key = "last point";
		else key = "first point";
		std::map<std::string,void*>::iterator it = msh->hook.find(key);
		if ( it != msh->hook.end() )     // found
			p = (Cell*) it->second;                               }
	if ( p == NULL )
	// p is still NULL, which means that the mesh has no "first point" attribute
	// (or "last point") thus we need to search for it
	// if the chain is not open (if it's a loop) this search may last forever !
	{	p = *(msh->cells[0]->begin());
	 	Cell *s = omsh->cell_in_front_of ( p, tag::may_not_exist );
	 	while ( s )
		{	p = & ( s->tip() );
	 		s = omsh->cell_in_front_of ( p, tag::may_not_exist );      }
		msh->hook[key] = (void*) p;                                        }
	Mesh::reset_iter_around_rev ( it, p );                                   }


void Mesh::reset_iter_around ( CellIterator *it, Cell * p )  // static
	
// if p is null, initialize the iterator at the beginning of the list of segments
{	Mesh::data_for_iter_contour *data =
		(Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->cells[msh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	Mesh *omsh = data->oriented_base;
	if ( p == NULL ) p = *(msh->cells[0]->begin());
  else if ( p->dim == 1 )
	{	Cell & s = *p;
		assert ( & ( omsh->orient_cell(s) ) == &s );
		p = s.base().hidden_reverse;
		assert ( p != NULL );                         }
	assert ( p->dim == 0 );
	assert ( p->is_positive() );
	data->current_vertex = p;
	data->current_seg = data->first_seg
	                  = omsh->cell_in_front_of ( p, tag::may_not_exist );
	data->first_try = ( data->first_seg != NULL );
	// data->first_try might be false if we come from reset_iter_along
	data->seg_got_out = false;                                                }


void Mesh::reset_iter_around_rev ( CellIterator *it, Cell * p )  // static
	
// if p is null, initialize the iterator at the beginning of the list of segments
{	Mesh::data_for_iter_contour *data =
		(Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	assert ( msh->cells[msh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	Mesh *omsh = data->oriented_base;
	if ( p == NULL ) p = *(msh->cells[0]->begin());
	else if ( p->dim == 1 )
	{	Cell & s = *p;
		assert ( & ( omsh->orient_cell(s) ) == &s );
		p = & ( s.tip() );                       }
	assert ( p->dim == 0 );
	data->current_vertex = p;
	data->current_seg = data->first_seg
	                  = omsh->cell_behind ( p, tag::may_not_exist );
	data->first_try = ( data->first_seg != NULL );
	// data->first_try might be false if we come from reset_iter_along
	data->seg_got_out = false;                                             }


void Mesh::advance_iter_contour ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	if ( data->seg_got_out ) // only meaningful for "vertices along"
	{	data->seg_got_out = false; return;  }
	data->first_try = false;
	Mesh *omsh = data->oriented_base;
	data->current_vertex = & ( data->current_seg->tip() );
	data->current_seg = omsh->cell_in_front_of
		( data->current_vertex, tag::may_not_exist );
	if ( data->current_seg == NULL ) data->seg_got_out = 1;                        }


void Mesh::advance_iter_contour_rev ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	if ( data->seg_got_out ) // only meaningful for "vertices along"
	{	data->seg_got_out = false; return;  }
	data->first_try = false;
	Mesh *omsh = data->oriented_base;
	data->current_vertex = data->current_seg->base().hidden_reverse;
	data->current_seg = omsh->cell_behind
		( data->current_vertex, tag::may_not_exist );
	if ( data->current_seg == NULL ) data->seg_got_out = 1;                          }


bool Mesh::valid_iter_around ( CellIterator *it )  // static

// the case when 'this' is a void mesh is not contemplated 

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	if ( data->first_try ) return 1;
	return ( data->current_seg != data->first_seg );                             }


bool Mesh::valid_iter_segs_along ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	return ( data->current_seg != NULL );                                        }


bool Mesh::valid_iter_vertices_along ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	if ( data->current_seg == NULL )  return ( data->seg_got_out );
	return 1;                                                                     }


Cell *Mesh::deref_iter_segs_contour ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	return data->current_seg;                                                     }


Cell *Mesh::deref_iter_vertices_contour ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	return data->current_vertex;                                                  }


void Mesh::reset_iter_given_dim ( CellIterator *it, Cell *p )  // static

{	assert ( p == NULL );
	Mesh::data_for_iter_given_dim *data =
		(Mesh::data_for_iter_given_dim*) it->data;
	Mesh *msh = (Mesh*) it->base;
	data->cell_list = msh->cells [data->dim];
	data->iter = data->cell_list->begin();        }


void Mesh::advance_iter_given_dim (CellIterator *it)  // static

{	Mesh::data_for_iter_given_dim *data =
		(Mesh::data_for_iter_given_dim*) it->data;
	data->iter++;                                 }


bool Mesh::valid_iter_given_dim (CellIterator *it)  // static

{	Mesh::data_for_iter_given_dim *data =
		(Mesh::data_for_iter_given_dim*) it->data;
	return ( data->iter != data->cell_list->end() ); }


Cell *Mesh::deref_iter_given_dim (CellIterator* it)  // static

{	Mesh::data_for_iter_given_dim *data =
		(Mesh::data_for_iter_given_dim*) it->data;
	Cell *cll = *data->iter;
	assert ( cll->is_positive() );
	return cll;                                    }

void Mesh::reset_iter_along_dual ( CellIterator *it, Cell *seg )  // static

// if 'seg' is NULL, initialize the iterator at the beginning of the chain

{	Mesh::data_for_iter_contour *data =
		(Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	Mesh *ambient_mesh = data->oriented_base;
	assert ( ambient_mesh->cells[ambient_mesh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	if ( seg == NULL )
	// we need to search for a starting point
	// if the cell 'data->center' is interior to ambient_mesh
	// (i.e. not on the boundary), this search may last forever !
	{	Cell *tri = ambient_mesh->get_tri_for_iter_dual ( data->center );
		// a triangle having 'data->center' as vertex, within ambient_mesh
		// what if 'data->center' is a negative cell ?
		Cell * P = data->center;
	 	while ( true )
	 	{	seg = tri->boundary().cell_in_front_of ( data->center, tag::may_not_exist );
	 		assert ( seg != NULL );
			Cell * rev_seg = seg->hidden_reverse;
			// 'rev_seg' (if not NULL) points towards 'data->center'
			// if seg has no reverse, we know it belongs to the boundary of ambient_mesh
			// if seg has a reverse, we cannot be sure - we must check further
			// because, for instance, there may be another, adjacent, mesh
			if ( rev_seg == NULL ) break;
			Cell * next_tri = ambient_mesh->cell_behind ( rev_seg, tag::may_not_exist );
			if ( next_tri == NULL ) break;  // 'seg' is on the boundary of ambient_mesh
			tri = next_tri;                                                                }
		seg = tri->boundary().cell_behind ( data->center, tag::may_not_exist );
		assert ( seg != NULL );                                                             }
	Mesh::reset_iter_around_dual ( it, seg );
}	// end of Mesh::reset_iter_along_dual


void Mesh::reset_iter_along_dual_rev ( CellIterator *it, Cell *seg )  // static

// if 'seg' is NULL, initialize the iterator at the beginning of the chain

{	Mesh::data_for_iter_contour *data =
		(Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	Mesh *ambient_mesh = data->oriented_base;
	assert ( ambient_mesh->cells[ambient_mesh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	if ( seg == NULL )
	// we need to search for a starting point
	// if the cell 'data->center' is interior to ambient_mesh
	// (i.e. not on the boundary), this search may last forever !
	{	Cell *tri = ambient_mesh->get_tri_for_iter_dual ( data->center );
		// a triangle having 'data->center' as vertex, within ambient_mesh
		// what if 'data->center' is a negative cell ?
	 	while ( true )
	 	{	seg = tri->boundary().cell_behind ( data->center, tag::may_not_exist );
	 		assert ( seg != NULL );
			// 'seg' points towards 'data->center'
			Cell * rev_seg = seg->hidden_reverse;
			// if seg has no reverse, we know it belongs to the boundary of ambient_mesh
			// if seg has a reverse, we cannot be sure - we must check further
			// because, for instance, there may be another, adjacent, mesh
			if ( rev_seg == NULL ) break;
			Cell * next_tri = ambient_mesh->cell_behind ( rev_seg, tag::may_not_exist );
			if ( next_tri == NULL ) break;  // 'seg' is on the boundary of ambient_mesh
			tri = next_tri;                                                               }
		seg = tri->boundary().cell_behind ( data->center, tag::may_not_exist );
		assert ( seg != NULL );                                                            }
	Mesh::reset_iter_around_dual ( it, seg );
}	// end of Mesh::reset_iter_along_dual_rev


void Mesh::reset_iter_around_dual ( CellIterator *it, Cell * seg )  // static

// if p is null, initialize the iterator at the beginning of the list of segments

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	Mesh *msh = (Mesh*) it->base;
	Mesh *ambient_mesh = data->oriented_base;
	assert ( ambient_mesh->cells[ambient_mesh->dim]->size() > 0 );
	// the case of a void mesh is not contemplated
	// if we arrived here from Mesh::reset_iter_along_dual, 'seg' is not NULL
	// but it could be a triangle
	if ( seg == NULL )
	{	Cell *tri = ambient_mesh->get_tri_for_iter_dual ( data->center );
		// a triangle having 'data->center' as vertex, within ambient_mesh
		// what if 'data->center' is a negative cell ?
		seg = tri->boundary().cell_behind ( data->center, tag::may_not_exist );  }
  else if ( seg->dim == data->center->dim + 2 )
	{	Cell * & tri = seg;
		assert ( & ( ambient_mesh->orient_cell(*tri) ) == tri );
		seg = tri->boundary().cell_behind ( data->center, tag::may_not_exist );  }
	assert ( seg != NULL );
	assert ( seg->dim == data->center->dim + 1 );
	assert ( & ( seg->boundary().orient_cell(*(data->center)) ) == data->center );
	data->current_vertex = seg;  // segment or higher; names misleading
	data->current_seg = data->first_seg  // triangle or higher; names misleading
	       = ambient_mesh->cell_behind ( seg, tag::may_not_exist );
	assert ( data->current_seg != NULL );
	data->first_try = true;
} // end of Mesh::reset_iter_around_dual


void Mesh::advance_iter_dual ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	Mesh *ambient_mesh = data->oriented_base;
	data->first_try = false;
	// names are misleading :
	// data->current_seg is not a segment now, it's a triangle or even higher dim
	// data->current_vertex is a face of data->current_seg
	// and data->center is a face of data->current_vertex (so, it's a subface)
	data->current_seg = ambient_mesh->cell_in_front_of
		( data->current_vertex, tag::may_not_exist );
	if ( data->current_seg == NULL ) return;
	//	data->seg_got_out = 1;
	data->current_vertex = data->current_seg->boundary().cell_behind
		( data->center, tag::may_not_exist );
	assert ( data->current_vertex != NULL );	                                }


void Mesh::advance_iter_dual_rev ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	Mesh *ambient_mesh = data->oriented_base;
	data->first_try = false;
	// names are misleading :
	// data->current_seg is not a segment now, it's a triangle or even higher dim
	// data->current_vertex is a face of data->current_seg
	// and data->center is a face of data->current_vertex (so, it's a subface)
	Cell * seg = data->current_seg->boundary().cell_in_front_of
		( data->center, tag::may_not_exist );
	assert ( seg != NULL );
	Cell * rev_seg = seg->hidden_reverse;
	if ( rev_seg == NULL )  // we're stuck
	{	data->current_seg = NULL;  return;  } 
	data->current_vertex = rev_seg;
	data->current_seg = ambient_mesh->cell_behind ( rev_seg, tag::may_not_exist );  }
	// if ( data->current_seg == NULL ) data->seg_got_out = 1;


bool Mesh::valid_iter_along_dual ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	return ( data->current_seg != NULL );                                         }


bool Mesh::valid_iter_around_dual ( CellIterator *it )  // static

{	Mesh::data_for_iter_contour *data = (Mesh::data_for_iter_contour*) it->data;
	if ( data->first_try ) return true;
	return ( data->current_seg != data->first_seg );                                 }


////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//      HHHHHH  H     H   H    H     HHH   HHHHHHH  H     HHH     H    H    HH        //
//      H       H     H   HH   H    H   H     H     H    H   H    HH   H   H  H       //
//      H       H     H   HH   H   H          H     H   H     H   HH   H   H          //
//      HHHH    H     H   H H  H   H          H     H   H     H   H H  H    H         //
//      H       H     H   H H  H   H          H     H   H     H   H H  H     H        //
//      H       H     H   H  H H   H          H     H   H     H   H  H H      H       //
//      H        H   H    H   HH    H   H     H     H    H   H    H   HH   H  H       //
//      H         HHH     H    H     HHH      H     H     HHH     H    H    HH        //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////


// factory functions :

FunctionOnMesh::baseFunction * FunctionOnMesh::sum
	( std::list < FunctionOnMesh::baseFunction * > l )
{	double c = 0.;  // we gather all constants here
	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	      it = l.begin(); it != l.end(); )
	{	if ( *it == FunctionOnMesh::zero.components[0] )
		{	it = l.erase (it);
			continue;           }
		if ( (*it)->type == "constant" )
		{	FunctionOnMesh::Constant * co = (FunctionOnMesh::Constant*) *it;
			c += co->c; 
			it = l.erase (it);
			continue;           }
		if ( (*it)->type == "sum" )
		{	FunctionOnMesh::Sum * s = (FunctionOnMesh::Sum*) *it;
			for (std::list<FunctionOnMesh::baseFunction*>::iterator
			     itt = s->terms.begin();
				itt != s->terms.end(); itt++ )
				l.push_back ( *itt ); 
			it = l.erase (it);
			continue;                                                  }
		it++;                                                                  }
	if ( c != 0. ) l.push_front ( new FunctionOnMesh::Constant(c) );
	if ( l.size() == 0 ) return FunctionOnMesh::zero.components[0];
	if ( l.size() == 1 ) return *(l.begin());
	return (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Sum(l);             }

FunctionOnMesh::baseFunction * FunctionOnMesh::product
	( std::list < FunctionOnMesh::baseFunction * > l )
{	double c = 1.;  // we gather all constants here
	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	      it = l.begin(); it != l.end(); )
	{	if ( *it == FunctionOnMesh::zero.components[0] )
			return FunctionOnMesh::zero.components[0];
		if ( *it == FunctionOnMesh::one.components[0] )
		{	it = l.erase (it);
			continue;           }
		if ( (*it)->type == "constant" )
		{	FunctionOnMesh::Constant * co = (FunctionOnMesh::Constant*) *it;
			c *= co->c; 
			it = l.erase (it);
			continue;                               }
		if ( (*it)->type == "product" )
		{	FunctionOnMesh::Product * p = (FunctionOnMesh::Product*) *it;
			for ( std::list<FunctionOnMesh::baseFunction*>::iterator
			      itt = p->factors.begin();	itt != p->factors.end(); itt++ )
				l.push_back ( *itt ); 
			it = l.erase (it);
			continue;                                                           }
		it++;                                                                    }
	if ( c != 1. ) l.push_front ( new FunctionOnMesh::Constant(c) );
	if ( l.size() == 0 ) return FunctionOnMesh::one.components[0];
	if ( l.size() == 1 ) return *(l.begin());
	return (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Product(l);         }

FunctionOnMesh::Function & FunctionOnMesh::from_field
	( NumericField & field, const std::string type )
{	assert ( type == "Lagrange degree one");
	size_t dim = field.size();
	if ( dim == 1 )
	{	assert ( field.type == "one-dim" );
		OneDimField & f = * ( (OneDimField*) &field );
		assert ( f.lives_on_cells_of_dim == 0 );        }
	else
	{	assert ( ( field.type == "wrapped block" ) ||
		         ( field.type == "wrapped composed" ) );
		MultiDimField * fi = (MultiDimField*) &field;
		if ( field.type == "wrapped block" )
		{	BlockField & f = * ( (BlockField*) fi->base );
			assert ( f.lives_on_cells_of_dim == 0 );      }
		else
		{	ComposedField & f = * ( (ComposedField*) fi->base );
			std::vector<BlockFieldBase*>::iterator it = f.fields.begin();
			for ( ; it != f.fields.end(); it++ )
				assert ( (*it)->lives_on_cells_of_dim == 0 );           }          }
	FunctionOnMesh::Function & f = * ( new FunctionOnMesh::Function ( dim ) );
	for ( size_t i = 0; i < dim; i++ )
		f.components[i] = new FunctionOnMesh::baseFunction;
	// we provide a NULL mesh, which means "for all cells, in any mesh" :
	AnyCell & cll = * ( new AnyCell ( (Mesh*) NULL, 0 ) );
	for ( size_t i = 0; i < dim; i++ ) f[i] (cll) == field[i] (cll);
	return f;                                                                    }

double FunctionOnMesh::baseFunction::operator() ( Cell & cll )
// we search between all relations involving 'this' variable
// if we find a relation, great
{	double result; short int found = 0;
	for ( std::list<RelationFunc*>::iterator it = this->relations.begin();
	      it != this->relations.end(); it++ )
  {	RelationFunc * rel = *it;
		assert ( rel->func == this );
		// check if the cell 'cll' is valid
		// ------------------------------------------------------------------------------------
		// we should take into account the case when 'cll' is part of ac->env_cell->boundary()
		// or maybe just a point inside ac->env_cell
		// ------------------------------------------------------------------------------------
		AnyCell *ac = rel->cell;
		if ( ac->env_cell )
		{	if ( ac->env_cell != &cll ) continue;  }
		else if ( ac->env_mesh )
		{	std::map < Mesh *, Cell::field_to_meshes > * cm =
				cll.meshes[ac->env_mesh->dim];
			if ( cm->find(ac->env_mesh) == cm->end() ) continue;  }
		// if both ac->env_cell and ac->env_mesh are NULL,
		// then 'ac' is valid in any mesh so no verification is necessary
		std::pair <bool,double> p = rel->operator()(cll);
		if ( p.first )  { found++; result = p.second;  }           }
	assert ( found == 1 ); return result;                                  }

double FunctionOnMesh::Function::operator() ( Cell & cll )
{	assert ( this->components.size() == 1 );
	FunctionOnMesh::baseFunction * f = *(this->components.begin());
	double result; short int found = 0;
	return (*f) (cll);                                               }

std::string FunctionOnMesh::baseFunction::repr ( std::string from )
{	return name;  }

std::string FunctionOnMesh::Function::repr ( std::string from )
{	assert ( components.size() == 1 );
	return components[0]->repr ( from );	}

std::string FunctionOnMesh::Constant::repr ( std::string from )
{	std::stringstream ss;
	ss << c; 
	std::string s = ss.str();
	if ( c < 0. ) s = "(" + s + ")";
	return s;                         }

std::string FunctionOnMesh::Sum::repr ( std::string from )
	{	std::list<FunctionOnMesh::baseFunction*>::iterator it = terms.begin();
		assert ( it != terms.end() );
		//std::cout << "-sum-repr- " << this << ", " << terms.size() << " items" << std::endl << flush;
		//std::cout << "-sum-repr- first " << (*it)->type << std::endl << flush;
		std::string s = (*it)->repr ("sum");
		//std::cout << "-sum-repr- "  << s << std::endl << flush;
		for ( it++; it != terms.end(); it++ )
		{	//std::cout << "-sum-repr- " << k <<"th " << (*it)->type <<flush;
			//std::cout << "-sum-repr- " << (*it)->repr() << std::endl << flush;
			s = s + '+' + (*it)->repr ("sum");  }
		if ( ( from == "power" ) || ( from == "product" ) ) s = "(" + s + ")";
		return s;                                                               }

std::string FunctionOnMesh::Product::repr ( std::string from )
	{	std::list<FunctionOnMesh::baseFunction*>::iterator it = factors.begin();
		assert ( it != factors.end() );
		//std::cout << "-prod-repr- " << this << ", " << factors.size() << " items" << std::endl << flush;
		//std::cout << "-prod-repr- first " << (*it)->type << std::endl << flush;
		std::string s = (*it)->repr ("product");
		//std::cout << "-prod-repr- " << s << std::endl << flush;
		for ( it++; it != factors.end(); it++ )
		{	//std::cout << "-prod-repr- " << k << "th " << (*it)->type << std::endl << flush;
			s = s + '*' + (*it)->repr ("product"); }
		if ( from == "power" ) s = "(" + s + ")";
		return s;                                                    }

std::string FunctionOnMesh::Fraction::repr ( std::string from )
{		std::string s = up->repr ("product");
		s = s + '/' + down->repr ("power");
		if ( from == "power" ) s = "(" + s + ")";
		return s;                                                    }

std::string FunctionOnMesh::Power::repr ( std::string from )
	{	std::string s = base->repr ("power") + "^";
		std::stringstream ss;
		ss << exponent;
		return s+ss.str();                      }

std::string FunctionOnMesh::Composition::repr ( std::string from )
{	return "(" + func->repr ( "composition" ) + ")oPhi";  }		 


void FunctionOnMesh::Mapping::build_jacobian ()

// The method "build_jacobian" computes symbolically the partial derivatives
// and stores them in "jacobian". If this matrix is square, it also computes
// the inverse matrix (as symbolic functions) and stores it in "inv_jacobian".

{	std::set < FunctionOnMesh::baseFunction* > vars_of_this;
	this->seek_vars ( vars_of_this );
	std::set<FunctionOnMesh::baseFunction*>::iterator it_var;
	for ( it_var = vars_of_this.begin(); it_var != vars_of_this.end(); it_var++ )
	{	FunctionOnMesh::baseFunction * variable = *it_var;
		std::vector<FunctionOnMesh::baseFunction*>::iterator it_comp;
		for ( it_comp = this->components.begin();
					it_comp != this->components.end(); it_comp++ )
		{	FunctionOnMesh::baseFunction * component = *it_comp;
			this->jacobian[component][variable] = component->deriv(variable);  }  }
	// we now want to compute the determinant of the jacobian and the inverse matrix
	// We distinguish two cases. First, when the jacobian is square (we have the
	// same number of components as variables). Then the determinant is computed
	// in the usual way and might be negative. Second, we may have more components
	// than variables, for instance if we work on a surface in 3D. Then there is no
	// determinant in the usual sense but, if there are only two variables, we can
	// compute the area-transformation coefficient as the square root of the sum of
	// (x,s y,t - x,t y,s)^2 (a sort of exterior product).
	// The case of more than two variables remains open.
	assert ( vars_of_this.size() == 2 );
	it_var = vars_of_this.begin();
	FunctionOnMesh::baseFunction * s = *it_var;
	it_var++;
	FunctionOnMesh::baseFunction * t = *it_var;
	if ( this->components.size() == 2 ) // first case
	{	std::vector<FunctionOnMesh::baseFunction*>::iterator it_comp =
			this->components.begin();
		FunctionOnMesh::baseFunction * x = *it_comp;
		it_comp++;
		FunctionOnMesh::baseFunction * y = *it_comp;
		FunctionOnMesh::baseFunction * dxds = this->jacobian[x][s];
		FunctionOnMesh::baseFunction * dxdt = this->jacobian[x][t];
		FunctionOnMesh::baseFunction * dyds = this->jacobian[y][s];
		FunctionOnMesh::baseFunction * dydt = this->jacobian[y][t];
		FunctionOnMesh::baseFunction * p1 = FunctionOnMesh::product ( dxds, dydt );
		FunctionOnMesh::baseFunction * p2 = FunctionOnMesh::product
			     ( FunctionOnMesh::minus_one.components[0], dxdt, dyds );
		this->det = FunctionOnMesh::sum ( p1, p2 );
		// matrix 'this->inv_jacobian' is stored oddly as a double map (column first)
		//  d s / d x :
		this->inv_jacobian[x][s] = FunctionOnMesh::fraction ( dydt, this->det );
		FunctionOnMesh::baseFunction * temp = FunctionOnMesh::product
			( FunctionOnMesh::minus_one.components[0], dxdt );
		//  d s / d y :
		this->inv_jacobian[y][s] = FunctionOnMesh::fraction ( temp, this->det );
		temp = FunctionOnMesh::product
			( FunctionOnMesh::minus_one.components[0], dyds );
		//  d t / d x :
		this->inv_jacobian[x][t] = FunctionOnMesh::fraction ( temp, this->det );
		//  d t / d y :
		this->inv_jacobian[y][t] = FunctionOnMesh::fraction ( dxds, this->det );     }
	else // second case
	// compute the area-transformation coefficient as the square root of the sum of
	// (x,s y,t - x,t y,s)^2 (a sort of exterior product).
	{	std::list<FunctionOnMesh::baseFunction*> list_of_sq;
		std::vector<FunctionOnMesh::baseFunction*>::iterator it1, it2;
		for ( it1 = this->components.begin(); it1 != this->components.end(); it1++ )
		{	FunctionOnMesh::baseFunction * x = *it1;
			FunctionOnMesh::baseFunction * dxds = this->jacobian[x][s];
			FunctionOnMesh::baseFunction * dxdt = this->jacobian[x][t];
			for ( it2 = it1, it2++; it2 != this->components.end(); it2++ )
			{	FunctionOnMesh::baseFunction * y = *it2;
				FunctionOnMesh::baseFunction * dyds = this->jacobian[y][s];
				FunctionOnMesh::baseFunction * dydt = this->jacobian[y][t];
				FunctionOnMesh::baseFunction * p1 = FunctionOnMesh::product ( dxds, dydt );
				FunctionOnMesh::baseFunction * p2 = FunctionOnMesh::product
					( FunctionOnMesh::minus_one.components[0], dxdt, dyds );
				FunctionOnMesh::baseFunction * dif = FunctionOnMesh::sum ( p1, p2 );
				FunctionOnMesh::baseFunction * dif2 = FunctionOnMesh::power ( dif, 2. );
				list_of_sq.push_back ( dif2 );                                              }  }
		FunctionOnMesh::baseFunction * sum_of = FunctionOnMesh::sum ( list_of_sq );
		this->det = FunctionOnMesh::power ( sum_of, 0.5 );                                   }
				
}

//////////////////////////////////////////////////////////////////////////////////
////////////////   variational problems   ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void VariationalProblem::discretize ( )

{	FunctionOnMesh::combinIntegrals * left = this->equation.lhs;
	FunctionOnMesh::combinIntegrals * right = this->equation.rhs;

	// discretize the left hand side :
	std::list < std::pair<double,FunctionOnMesh::Integral*> > ::iterator it_integr;
	for ( it_integr = left->terms.begin(); it_integr != left->terms.end(); it_integr++ )
	{	FunctionOnMesh::Integral * integr = it_integr->second;
		Mesh * domain = integr->domain;
		CellIterator it = domain->iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( it.reset(); it.in_range(); it++ )
		{	Cell & cll = *it;
			this->finite_element->dock_on ( cll );
			this->add_to_matrix ( integr );         }                         }
	
	// discretize the right hand side :
//	for ( it_integr = right->terms.begin(); it_integr != right->terms.end(); it_integr++ )
//	{	FunctionOnMesh::Integral * integr = it_integr->second;
//		Mesh * domain = integr->domain;
//		CellIterator it = domain->iter_over ( tag::cells, tag::of_max_dim, tag::oriented );

} // end of VariationalProblem::discretize


void VariationalProblem::add_to_matrix ( FunctionOnMesh::Integral * integr )

// add contribution of current cell to the global matrix

// podemos ir ainda mais a fundo e acrescentar, para cada ponto de Gauss,
// a contribuicao de todos os integrais para a matriz e vector globais
// o codigo talvez fique mais rapido assim ...

{	FiniteElement * fe = this->finite_element;
	Cell * cll = FiniteElement::current_cell;
	Integrator & quadrature = * ( this->finite_element->integrator.at(cll->dim) );

	// if a master element is used, we need to replace
	// each space_coord with the corresponding composed_space_coord
	// (if we use no master element, do nothing)
	integr = fe->replace_space_coords ( integr );

	// then we replace the unknown by a local_base function
	// and the test by another local_base function
	FunctionOnMesh::baseFunction * u = * ( this->unknown.begin() );
	FunctionOnMesh::baseFunction * t = * ( this->test.begin() );
	std::list < FunctionOnMesh::baseFunction* > & psi = fe->local_base;
	std::list<FunctionOnMesh::baseFunction*>::iterator it1, it2;
	for ( it1 = psi.begin(); it1 != psi.end(); it1++ )
	{	// replace the unknown function with one psi
		FunctionOnMesh::Integral * integr_1 = integr->replace ( u, *it1 );
		for ( it2 = psi.begin(); it2 != psi.end(); it2++ )
		{	// replace the test function with another psi
			FunctionOnMesh::Integral * integr_2 = integr_1->replace ( t, *it2 );
			// usually, there are delayed derivatives; compute them
			integr_2 = integr_2->compute_derivatives();
			std::cout << quadrature ( * ( integr_2->integrand ) ) << " ";
			// add_to_matrix_at_degrees_of_freedom;
		} std::cout << std::endl; }  }

//////////////////////////////////////////////////////////////////////////////////
///////////////    finite elements, integrators   ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


FunctionOnMesh::baseFunction * FiniteElement::withMaster::replace_space_coords
( FunctionOnMesh::baseFunction * f )

// For finite elements which make use of a master element, we have two kinds of
// space coordinates. We have the "original" ones, say, x and y, and then we have
// something like x_comp = x(A)*psi_A + x(B)*psi_B + ..., y_comp = y(A)*psi_A + ...
// From the mathematical point of view, x and x_comp are one and the same function,
// because they take always the same values.
// However, C++ sees them as distinct FunctionOnMesh::baseFunction objects.
// The variational formulation has been defined before the finite element,
// in particular before the docking on a certain cell,
// and the integrals (and derivatives) appearing there are in terms of x and y.
// But then we want to integrate them in x_comp and y_comp through the transformation
// of coordinates from the physical element to the master element, so we need
// to replace, in each integral, x by x_comp and y by y_comp.

// For elements which use no master, this operation is not necessary, see
// the method FiniteElement::replace_space_coords in Mesh.h

{	size_t n = this->composed_space_coords.size();
	assert ( this->space_coords.size() == n );
	std::list<FunctionOnMesh::baseFunction*>::iterator
		it_space = this->space_coords.begin(),
		it_comp = this->composed_space_coords.begin();
	for ( size_t i = 0; i < n; i++ )
	{	assert ( it_space != this->space_coords.end() );
		assert ( it_comp != this->composed_space_coords.end() );
		f = f->replace_var ( * it_space, * it_comp );
		it_space++; it_comp++;                                     }
	assert ( it_space == this->space_coords.end() );
	assert ( it_comp == this->composed_space_coords.end() );
	return f;                                                        }


double hidden::GaussQuadrature
( Integrator * integr, FunctionOnMesh::baseFunction & f )

{	hidden::GaussIntegrator * gauss_int = (hidden::GaussIntegrator*) integr;
	std::list<FunctionOnMesh::baseFunction*>::iterator it = gauss_int->variables.begin();
	assert ( it != gauss_int->variables.end() );
	FunctionOnMesh::baseFunction & xi  = * ( *it );
	it++; assert ( it != gauss_int->variables.end() );
	FunctionOnMesh::baseFunction & eta = * ( *it );
	it++; assert ( it == gauss_int->variables.end() );
	FunctionOnMesh::baseFunction & jacobian = * ( gauss_int->jacobian );
	double result = 0.;
	size_t npg = gauss_int->gauss_point_coord.size();
	assert ( npg == gauss_int->gauss_weight.size() );
	for ( size_t i = 0; i < npg; i++)
	{	Cell & gauss_point = Cell::point();
		// this is not a vertex, but a point inside the cell
		// should we use a distinct class ?
		FunctionOnMesh::prescribe_on (gauss_point);
		xi  == gauss_int->gauss_point_coord[i][0];
		eta == gauss_int->gauss_point_coord[i][1];
		result += f ( gauss_point )
			* gauss_int->gauss_weight[i] * jacobian ( gauss_point );  }
	return result;                                                              }


void FiniteElement::set_integrator
( std::string s1, short int d, std::string s2, short int n, std::string s3 )

{	assert ( s1 == "on cells of dimension");
	assert ( s2 == "gauss");
	assert ( s3 == "nodes");
	assert ( n == 9 );
	assert ( d == 2 );
	this->integrator[d] = & Integrator::gauss ("Q9", "set integrator");
}


//////////////////////////////////////////////////////////////////////////////////
////////////////   specific finite elements   ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


FiniteElement::Lagrange & FiniteElement::lagrange ( std::string s )

// this version with one argument is called before building the mesh
// reserves space for enumeration of future cells

{	FiniteElement::Lagrange * fe;
	if ( ( s == "Q1" ) || ( s == "Q1, no enum" ) )
		fe = new FiniteElement::Lagrange::Q1;
	else if ( ( s == "P1" ) || ( s == "P1, no enum" ) )
		assert (0);  // fe = new FiniteElement::Lagrange::P1;
	else if ( ( s == "Q2" ) || ( s == "Q2, no enum" ) )
		fe = new FiniteElement::Lagrange::Q2;
	else if ( ( s == "P2" ) || ( s == "P2, no enum" ) )
		assert (0);  // fe = new FiniteElement::Lagrange::P2;
	else assert ( ( s == "P1, no enum" ) || ( s == "Q1, no enum" ) );
	if ( ( s == "Q1" ) || ( s == "P1" ) )
	{	// we want to enumerate future vertices
		FiniteElement::Lagrange::degree_1 * f_e =
			( FiniteElement::Lagrange::degree_1 * ) fe;
		f_e->max_dof = 0;
		f_e->pointer_into_heap = Cell::size_t_heap_size[0];
		Cell::size_t_heap_size[0]++;
		Cell::init_cell[0]->push_back
			( & FiniteElement::Lagrange::degree_1::enumerate_new_vertex );
		Cell::data_for_init[0]->push_back ( (void*) fe );                }
	else if ( ( s == "Q2" ) || ( s == "P2" ) )
	{	// we want to enumerate future vertices and segments
		FiniteElement::Lagrange::degree_2 * f_e =
			( FiniteElement::Lagrange::degree_2 * ) fe;
		f_e->max_dof = 0;
		f_e->pointer_into_heap[0] = Cell::size_t_heap_size[0];
		Cell::size_t_heap_size[0]++;
		Cell::init_cell[0]->push_back
			( & FiniteElement::Lagrange::degree_2::enumerate_new_vertex );
		Cell::data_for_init[0]->push_back ( (void*) fe );
		f_e->pointer_into_heap[1] = Cell::size_t_heap_size[1];
		Cell::size_t_heap_size[1]++;
		Cell::init_cell[1]->push_back
			( & FiniteElement::Lagrange::degree_2::enumerate_new_segment );
		Cell::data_for_init[1]->push_back ( (void*) fe );
		if ( s == "Q2" ) // we also want to enumerate quadrilaterals
		{	f_e->pointer_into_heap[2] = Cell::size_t_heap_size[2];
			Cell::size_t_heap_size[2]++;
			Cell::init_cell[2]->push_back
				( & FiniteElement::Lagrange::Q2::enumerate_new_quadr );
			Cell::data_for_init[2]->push_back ( (void*) fe );                  }  }
	// if no enumeration is required, do nothing
	return * fe;                                                                     }


FiniteElement::Lagrange & FiniteElement::lagrange ( std::string s, Mesh & msh )

// this version with two arguments is called when the mesh already exists
// re-sizes the size_t heap of existing cells for enumeration

{	FiniteElement::Lagrange & fe = FiniteElement::lagrange (s);
	if ( ( s == "Q1" ) || ( s == "P1" ) )
	{	// apparently a mesh already exists, we enumerate its vertices
		FiniteElement::Lagrange::degree_1 * f_e =
			( FiniteElement::Lagrange::degree_1 * ) (&fe);
		CellIterator it = msh.iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
		for ( it.reset(); it.in_range(); it++ )
		{	Cell & p = *it;
			assert ( p.size_t_heap.size() <= f_e->pointer_into_heap );
			p.size_t_heap.resize(f_e->pointer_into_heap+1);
			p.size_t_heap[f_e->pointer_into_heap] = f_e->max_dof;
			f_e->max_dof++;                                              }                }
	else if ( ( s == "Q2" ) || ( s == "P2" ) )
	{	// apparently a mesh already exists, we enumerate its vertices and segments
		FiniteElement::Lagrange::degree_2 * f_e =
			( FiniteElement::Lagrange::degree_2 * ) (&fe);
		CellIterator it0 = msh.iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
		for ( it0.reset(); it0.in_range(); it0++ )
		{	Cell & p = *it0;
			assert ( p.size_t_heap.size() <= f_e->pointer_into_heap[0] );
			p.size_t_heap.resize(f_e->pointer_into_heap[0]+1);
			p.size_t_heap[f_e->pointer_into_heap[0]] = f_e->max_dof;
			f_e->max_dof++;                                                    }
		CellIterator it1 = msh.iter_over ( tag::cells, tag::of_dim, 1, tag::not_oriented );
		for ( it1.reset(); it1.in_range(); it1++ )
		{	Cell & s = *it1;
			assert ( s.size_t_heap.size() <= f_e->pointer_into_heap[1] );
			s.size_t_heap.resize(f_e->pointer_into_heap[1]+1);
			s.size_t_heap[f_e->pointer_into_heap[1]] = f_e->max_dof;
			f_e->max_dof++;                                                    }
		if ( ( s == "Q2" ) ) // we also enumerate quadrilaterals
		{	CellIterator it2 = msh.iter_over ( tag::cells, tag::of_dim, 2, tag::not_oriented );
			for ( it2.reset(); it2.in_range(); it2++ )
			{	Cell & q = *it2;
				assert ( q.size_t_heap.size() <= f_e->pointer_into_heap[2] );
				q.size_t_heap.resize(f_e->pointer_into_heap[2]+1);
				q.size_t_heap[f_e->pointer_into_heap[2]] = f_e->max_dof;
				f_e->max_dof++;                                                    }          }  }
	else assert ( ( s == "Q1, no enum" ) || ( s == "Q2, no enum" ) );
	return  fe;                                                                               }


void FiniteElement::Lagrange::Q1::dock_on ( Cell & cll, std::string s ) // virtual

{	assert ( cll.dim == 2 );
	assert ( s == "dofs only");
	CellIterator itt = cll.boundary().iter_over ( tag::segments, tag::around );
	itt.reset();
	Cell & AB = *itt; itt++; assert (itt.in_range());
	Cell & BC = *itt; itt++; assert (itt.in_range());
	Cell & CD = *itt; itt++; assert (itt.in_range());
	Cell & DA = *itt; itt++; assert (!(itt.in_range()));
	Cell & A = DA.tip(), & B = AB.tip(), & C = BC.tip(), & D = CD.tip();
	this->dofs.clear();
	this->dofs.push_back ( &A );
	this->dofs.push_back ( &B );
	this->dofs.push_back ( &C );
	this->dofs.push_back ( &D );
}
	
void FiniteElement::Lagrange::Q1::dock_on ( Cell & cll ) // virtual

{	assert ( cll.dim == 2 );
	FiniteElement::current_cell = &cll;
	CellIterator itt = cll.boundary().iter_over ( tag::segments, tag::around );
	itt.reset();
	Cell & AB = *itt; itt++; assert (itt.in_range());
	Cell & BC = *itt; itt++; assert (itt.in_range());
	Cell & CD = *itt; itt++; assert (itt.in_range());
	Cell & DA = *itt; itt++; assert (!(itt.in_range()));
	Cell & A = DA.tip(), & B = AB.tip(), & C = BC.tip(), & D = CD.tip();
	this->dofs.clear();
	this->dofs.push_back ( &A );
	this->dofs.push_back ( &B );
	this->dofs.push_back ( &C );
	this->dofs.push_back ( &D );

	FunctionOnMesh::Function & xi = this->new_coord(), & eta = this->new_coord();
	xi .components[0]->name = "xi";
	eta.components[0]->name = "eta";
	FunctionOnMesh::prescribe_on (DA);  xi  == -1.;
  FunctionOnMesh::prescribe_on (BC);  xi  ==  1.;
  FunctionOnMesh::prescribe_on (AB);  eta == -1.;
  FunctionOnMesh::prescribe_on (CD);  eta ==  1.;
	// oddly enough, the assertions above are not actually used in the code
	// they would be useful if we wanted to compute, for instance, x_comp(A)
	// which should produce the same value as x(A)
	// but we will never need x_comp(A)

	// master coordinates may have remained from a previous dock :
	this->master_coords.clear();
	this->master_coords.push_back ( xi .components[0] );
	this->master_coords.push_back ( eta.components[0] );
	this->integrator.at(cll.dim)->variables = this->master_coords;

	// FunctionOnMesh::Function & quarter = FunctionOnMesh::constant ( 0.25 );
	FunctionOnMesh::Function & psiA = (1.-xi)*(1.-eta), & psiB = (1.+xi)*(1.-eta),
		                       & psiC = (1.+xi)*(1.+eta), & psiD = (1.-xi)*(1.+eta);
	// each takes value 4 in the associated vertex
	
	//  std::cout << "dofs : " << A.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << B.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << C.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << D.size_t_heap.at(this->pointer_into_heap) << std::endl;

	// list<RelationFunc*>::iterator it_rel = xi.components[0]->relations.begin();
	// RelationFuncConst* r = (RelationFuncConst*) *it_rel;
	// std::cout << &DA << " " << r->cell->env_mesh << " " << r->cell->env_cell << std::endl;
	// std::cout << &A << std::endl;
	// std::cout << xi(DA) << std::endl;

	// we now declare x and y (and z) as functions of xi and eta
	FunctionOnMesh::Function * T = new FunctionOnMesh::Function(0);
	T->components.reserve ( this->space_coords.size() );
	// composed space coordinates may have remained from a previous dock :
	this->composed_space_coords.clear();
	std::list<FunctionOnMesh::baseFunction*>::iterator it;
	for ( it = this->space_coords.begin(); it != this->space_coords.end(); it++ )
	{	FunctionOnMesh::Function x ( * ( *it ) );
		FunctionOnMesh::Function & x_comp = x(A)*psiA + x(B)*psiB + x(C)*psiC + x(D)*psiD;
		x_comp.components[0]->name = "comp";
		// do we need both T and composed_space_coords ?
		T->components.insert ( T->components.end(), x_comp.components[0] );
		this->composed_space_coords.push_back ( x_comp.components[0] );                     }
		
	// x and x_comp are mathematically identical (take always the same values)
	// see the instructions 'prescribe_on' above
	// but from C++'s viewpoint they are different objects
	// it would be nice to be able to differentiate psi_A with respect to x
	// but C++ only allows to differentiate psi_A with respect to x_comp
	// it's up to the human to note that x and x_comp can be identified
	// thus, for each integral which is evaluated, we will need to replace
	// each space_coord with the corresponding composed_space_coord
	// (the method 'replace_space_coords' does this)

	// pack them together :
	FunctionOnMesh::Mapping & map_T = FunctionOnMesh::mapping ( *T );
	Integrator::Gauss * gauss_int = (Integrator::Gauss*) this->integrator.at(cll.dim);
	gauss_int->jacobian = map_T.det;
	if ( this->space_coords.size() == this->master_coords.size() )
	{	// local base functions may have remained from a previous dock :
		this->local_base.clear();
		this->local_base.push_back
			( psiA .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiB .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiC .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiD .compose_with_inverse_of ( map_T ) .components[0] ); }
		// each gives value 4 in the associated vertex
	// else (if we have more space coordinates than master coordinates) do nothing

} // end of FiniteElement::Lagrange::Q1::dock_on


void FiniteElement::Lagrange::Q2::dock_on ( Cell & cll, std::string s ) // virtual

{	assert ( cll.dim == 2 );
	assert ( s == "dofs only");
	CellIterator itt = cll.boundary().iter_over ( tag::segments, tag::around );
	itt.reset();
	Cell & AB = *itt; itt++; assert (itt.in_range());
	Cell & BC = *itt; itt++; assert (itt.in_range());
	Cell & CD = *itt; itt++; assert (itt.in_range());
	Cell & DA = *itt; itt++; assert (!(itt.in_range()));
	Cell & A = DA.tip(), & B = AB.tip(), & C = BC.tip(), & D = CD.tip();
	this->dofs.clear();
	this->dofs.push_back ( &A );
	this->dofs.push_back ( &B );
	this->dofs.push_back ( &C );
	this->dofs.push_back ( &D );
	this->dofs.push_back ( &AB );
	this->dofs.push_back ( &BC );
	this->dofs.push_back ( &CD );
	this->dofs.push_back ( &DA );
	this->dofs.push_back ( &cll );
}
	
void FiniteElement::Lagrange::Q2::dock_on ( Cell & cll ) // virtual

{	assert ( cll.dim == 2 );
	FiniteElement::current_cell = &cll;
	
	CellIterator itt = cll.boundary().iter_over ( tag::segments, tag::around );
	itt.reset();
	Cell & AB = *itt; itt++; assert (itt.in_range());
	Cell & BC = *itt; itt++; assert (itt.in_range());
	Cell & CD = *itt; itt++; assert (itt.in_range());
	Cell & DA = *itt; itt++; assert (!(itt.in_range()));
	Cell & A = DA.tip(), & B = AB.tip(), & C = BC.tip(), & D = CD.tip();
	this->dofs.clear();
	this->dofs.push_back ( &A );
	this->dofs.push_back ( &B );
	this->dofs.push_back ( &C );
	this->dofs.push_back ( &D );
	this->dofs.push_back ( &AB );
	this->dofs.push_back ( &BC );
	this->dofs.push_back ( &CD );
	this->dofs.push_back ( &DA );
	this->dofs.push_back ( &cll );

	FunctionOnMesh::Function & xi = this->new_coord(), & eta = this->new_coord();
	xi .components[0]->name = "xi";
	eta.components[0]->name = "eta";
	FunctionOnMesh::prescribe_on (DA);  xi  == -1.;
  FunctionOnMesh::prescribe_on (BC);  xi  ==  1.;
  FunctionOnMesh::prescribe_on (AB);  eta == -1.;
  FunctionOnMesh::prescribe_on (CD);  eta ==  1.;
	// oddly enough, the assertions above are not actually used in the code
	// they would be useful if we wanted to compute, for instance, x_comp(A)
	// which should produce the same value as x(A)
	// but we will never need x_comp(A)

	// master coordinates may have remained from a previous dock :
	this->master_coords.clear();
	this->master_coords.push_back ( xi .components[0] );
	this->master_coords.push_back ( eta.components[0] );
	this->integrator.at(cll.dim)->variables = this->master_coords;

	// FunctionOnMesh::Function & quarter = FunctionOnMesh::constant ( 0.25 );
	FunctionOnMesh::Function & psiA = (1.-xi)*(1.-eta), & psiB = (1.+xi)*(1.-eta),
		                       & psiC = (1.+xi)*(1.+eta), & psiD = (1.-xi)*(1.+eta);
	// each gives value 4 in the associated vertex
	
	//  std::cout << "dofs : " << A.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << B.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << C.size_t_heap.at(this->pointer_into_heap) << " "
	//  		 << D.size_t_heap.at(this->pointer_into_heap) << std::endl;

	// list<RelationFunc*>::iterator it_rel = xi.components[0]->relations.begin();
	// RelationFuncConst* r = (RelationFuncConst*) *it_rel;
	// std::cout << &DA << " " << r->cell->env_mesh << " " << r->cell->env_cell << std::endl;
	// std::cout << &A << std::endl;
	// std::cout << xi(DA) << std::endl;

	// we now declare x and y (and z) as functions of xi and eta
	FunctionOnMesh::Function * T = new FunctionOnMesh::Function(0);
	T->components.reserve ( this->space_coords.size() );
	// composed space coordinates may have remained from a previous dock :
	this->composed_space_coords.clear();
	std::list<FunctionOnMesh::baseFunction*>::iterator it;
	for ( it = this->space_coords.begin(); it != this->space_coords.end(); it++ )
	{	FunctionOnMesh::Function x ( * ( *it ) );
		FunctionOnMesh::Function & x_comp = x(A)*psiA + x(B)*psiB + x(C)*psiC + x(D)*psiD;
		x_comp.components[0]->name = "comp";
		// do we need both T and composed_space_coords ?
		T->components.insert ( T->components.end(), x_comp.components[0] );
		this->composed_space_coords.push_back ( x_comp.components[0] );                     }
		
	// x and x_comp are mathematically identical (take always the same values)
	// see the instructions 'prescribe_on' above
	// but from C++'s viewpoint they are different objects
	// it would be nice to be able to differentiate psi_A with respect to x
	// but C++ only allows to differentiate psi_A with respect to x_comp
	// it's up to the human to note that x and x_comp can be identified
	// thus, for each integral which is evaluated, we will need to replace
	// each space_coord with the corresponding composed_space_coord
	// (the method 'replace_space_coords' does this)

	// pack them together
	// the order in this->dofs should match exactly the order in this->local_base
	FunctionOnMesh::Mapping & map_T = FunctionOnMesh::mapping ( *T );
	Integrator::Gauss * gauss_int = (Integrator::Gauss*) this->integrator.at(cll.dim);
	gauss_int->jacobian = map_T.det;
	if ( this->space_coords.size() == this->master_coords.size() )
	{	// local base functions may have remained from a previous dock :
		this->local_base.clear();
		this->local_base.push_back
			( psiA .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiB .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiC .compose_with_inverse_of ( map_T ) .components[0] );
		this->local_base.push_back
			( psiD .compose_with_inverse_of ( map_T ) .components[0] );
		// each gives value 4 in the associated vertex
		this->local_base.push_back
			( (psiA*psiB) .compose_with_inverse_of(map_T) .components[0] );
		this->local_base.push_back
			( (psiB*psiC) .compose_with_inverse_of(map_T) .components[0] );
		this->local_base.push_back
			( (psiC*psiD) .compose_with_inverse_of(map_T) .components[0] );
		this->local_base.push_back
			( (psiD*psiA) .compose_with_inverse_of(map_T) .components[0] );
		// each gives value 4 in the center of the associated segment
		this->local_base.push_back
			( (psiA*psiB*psiC*psiD) .compose_with_inverse_of ( map_T ) .components[0] );
		// gives value 1 in the center of the quadrilateral
	}	// else (if we have more space coordinates than master coordinates) do nothing

} // end of FiniteElement::Lagrange::Q2::dock_on

