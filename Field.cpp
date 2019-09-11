// src/manifem/Field.cpp 2019.09.11

#include <cmath>
#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "assert.h"
#include "Mesh.h"


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//         HHHHHHH    H    HHHHHHH    H          HHHH        HHHH           //
//         H          H    H          H          H   HH     H    H          //
//         H          H    H          H          H     H    H               //
//         HHHHH      H    HHHHH      H          H     H      HHH           //
//         H          H    H          H          H     H         H          //
//         H          H    H          H          H     H          H         //
//         H          H    H          H          H   HH     H   HH          //
//         H          H    HHHHHHH    HHHHHHHH   HHHH        HHH            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


// Fields implement the idea of coordinates.
// They can also hold other numeric values like the solution of a PDE.

///////////////////////////////////////////////////////////////
//////////////     used by constructors    ////////////////////
///////////////////////////////////////////////////////////////

void MultiDimField::hidden_build ( BlockFieldBase *b )
{	this->base = (void*) b;
	this->restr = MultiDimField::restr_block;
	this->proj = MultiDimField::proj_block;
	this->join_this_block = MultiDimField::join_wrappedblock_block;
	this->join_this_composed = MultiDimField::join_wrappedblock_composed;
	this->join_this_multi = MultiDimField::join_wrappedblock_multi;
	this->join_block_this = MultiDimField::join_block_wrappedblock;
	this->join_composed_this = MultiDimField::join_composed_wrappedblock;
	this->join_multi_this = MultiDimField::join_multi_wrappedblock;
	this->interpolate2 = MultiDimField::interpolate_block_2;
	this->interpolate3 = MultiDimField::interpolate_block_3;
	this->interpolate4 = MultiDimField::interpolate_block_4;
	this->interpolate6 = MultiDimField::interpolate_block_6;                 }

void MultiDimField::hidden_build ( ComposedField *b )
{	this->base = (void*) b;
	this->restr = MultiDimField::restr_composed;
	this->proj = MultiDimField::proj_composed;
	this->join_this_block = MultiDimField::join_wrappedcomposed_block;
	this->join_this_composed = MultiDimField::join_wrappedcomposed_composed;
	this->join_this_multi = MultiDimField::join_wrappedcomposed_multi;
	this->join_block_this = MultiDimField::join_block_wrappedcomposed;
	this->join_composed_this = MultiDimField::join_composed_wrappedcomposed;
	this->join_multi_this = MultiDimField::join_multi_wrappedcomposed;
	this->interpolate2 = MultiDimField::interpolate_composed_2;
	this->interpolate3 = MultiDimField::interpolate_composed_3;
	this->interpolate4 = MultiDimField::interpolate_composed_4;
	this->interpolate6 = MultiDimField::interpolate_composed_6;                }


///////////////////////////////////////////////////////////////
///////////////    pseudo-constructors    /////////////////////
///////////////////////////////////////////////////////////////


OneDimField & NumericField::one_dim ( std::string s2, std::string type )
// examples of use :
//  auto &x = NumericField::one_dim ("lives on", "points");
//  auto &x = NumericField::one_dim ("lives on", "segments");
//  auto &x = NumericField::one_dim ("lives on", "2d cells");
//  auto &x = NumericField::one_dim ("lives on", "3d cells");

{	assert ( s2 == "lives on" );
	if ( type == "points" )
		return * ( new OneDimField (0) );
	if ( type == "segments" )
		return * ( new OneDimField (1) );
	if ( type == "2d cells" )
		return * ( new OneDimField (2) );
	if ( type == "3d cells" )
		return * ( new OneDimField (3) );    }


MultiDimField & NumericField::multi_dim ( std::string s1, size_t n, std::string s2, std::string type )
// examples of use :
//  auto &xy = NumericField::multi_dim ("size", 2, "lives on", "points");
//  auto &xy = NumericField::multi_dim ("size", 2, "lives on", "segments");
//  auto &xyz = NumericField::multi_dim ("size", 3, "lives on", "2d cells");
//  auto &xy = NumericField::multi_dim ("size", 2, "lives on", "3d cells");

{	assert ( s1 == "size" );
	assert ( s2 == "lives on" );
	if ( type == "points" )
		return * (new MultiDimField (n, 0) );
	if ( type == "segments" )
		return * (new MultiDimField (n, 1) );
	if ( type == "2d cells" )
		return * (new MultiDimField (n, 2) );
	if ( type == "3d cells" )
		return * (new MultiDimField (n, 3) );   }


///////////////////////////////////////////////////////////////
/////////////   size (number of components)   /////////////////
///////////////////////////////////////////////////////////////


// we forbid (at run time) instantiation of NumericField objects
size_t NumericField::size () const {  assert (0); };

// we forbid (at run time) instantiation of BlockFieldBase objects
size_t BlockFieldBase::size () const {  assert (0); };

size_t OneDimField::size() const {	return (size_t) 1;  }

size_t BlockField::size() const
{	return this->index_max_plus_one - this->index_min;  }

size_t ComposedField::size() const
{	size_t s = (size_t) 0;
	std::vector <BlockFieldBase*> :: const_iterator
		it = this->fields.begin(),
		it_e = this->fields.end();
	for ( ; it != it_e; it++ )  s += (*it)->size();
	return s;                                         }
	
size_t MultiDimField::size() const
{	return ( (NumericField*) this->base )->size();  }


///////////////////////////////////////////////////////////////
//////////////         operator ()        /////////////////////
///////////////////////////////////////////////////////////////
//////////////     coordinate of a cell   /////////////////////
///////////////////////////////////////////////////////////////


PointerVector<double> & BlockFieldBase::operator() (Cell &cll)

{	assert ( this->lives_on_cells_of_dim == cll.dim );
	assert ( cll.real_heap.size() >= this->index_max_plus_one );
	PointerVector<double> *res = new PointerVector<double> ( this->size() ) ;
	std::vector <double*> :: iterator it = res->pointers.begin();
	for ( size_t i = this->index_min; i < this->index_max_plus_one; i++, it++ )
		*it = &(cll.real_heap[i]);  // optimizar
	return *res;                                                                  }


PointerVector<double> & ComposedField::operator() (Cell &cll)

{	PointerVector<double> *res = new PointerVector<double> ( this->size() );
	std::vector <double*> :: iterator it_res = res->pointers.begin();
	std::vector <BlockFieldBase*> :: iterator itvf = this->fields.begin(),
		itvf_e = this->fields.end();
	for ( ; itvf != itvf_e; itvf++ )
	{	if ( (*itvf)->type == "one-dim" )
		{	OneDimField *f = (OneDimField*) *itvf;
			*it_res = &((*f)(cll));
			it_res++;                                }
		else
		{	assert ( (*itvf)->type == "block" );
			BlockField *f = (BlockField*) *itvf;
			std::vector <double*> &part = (*f)(cll).pointers;
			std::vector <double*> :: iterator it_part = part.begin(),
				it_part_e = part.end();
			for ( ; it_part != it_part_e; it_part++, it_res++ )
			{	assert ( it_res != res->pointers.end() );
				*it_res = *it_part;                             }  }  }
	assert ( it_res == res->pointers.end() );
	return *res;                                                                      }


///////////////////////////////////////////////////////////////
////////////     operator[] -- projection    //////////////////
///////////////////////////////////////////////////////////////

// used to get one component from a pair or a triplet

OneDimField & NumericField::operator[] ( size_t i )
{//	if ( this->type == "one-dim" )
  //	{	OneDimField *that = (OneDimField*) this;
  //		return that->operator[] (i);              }
  //	if ( this->type == "multi-dim" )
  //	{	MultiDimField *that = (MultiDimField*) this;
  //		return that->operator[] (i);                   }
	assert (0);                                               }

OneDimField & OneDimField::operator[] ( size_t i )
{	assert ( i == 0 );
	return *this;        }

OneDimField & MultiDimField::operator[] ( size_t i )
{	OneDimField &f = this->proj (this->base, i);
	f.is_part_of.insert ( (NumericField*) this );
	return f;                                     }

OneDimField & BlockField::operator[] ( size_t i )
{	assert ( i < this->size() );
	OneDimField *f = new OneDimField;
	f->lives_on_cells_of_dim = this->lives_on_cells_of_dim;
	f->index_min = this->index_min + i;
	f->index_max_plus_one = f->index_min + 1;
	f->is_part_of.insert (this);
	return *f;                                                }

OneDimField & ComposedField::operator[] ( size_t i )
{	assert ( i < this->size() );
	size_t s = (size_t) 0, old_s;
	std::vector <BlockFieldBase*> :: iterator it = this->fields.begin(),
		it_e = this->fields.end();
	for ( ; it != it_e; it++ )
	{	old_s = s; s += (*it)->size();
		if ( i < s ) break;                            }
	OneDimField &f = (*it)->operator[] ( i - old_s );
	assert ( f.is_part_of.find(this) != f.is_part_of.end() );
	return f;                                                                  }


///////////////////////////////////////////////////////////////
///////////////////     operator &&     ///////////////////////
///////////////////////////////////////////////////////////////

// used to "concatenate" two fields

// to do : update 'is_part_of' flag

MultiDimField & NumericField::operator&& ( NumericField &f )
{//	if ( this->type == "one-dim" )
  //	{	OneDimField *that = (OneDimField*) this;
  //		return that->operator&& (f);              }
  //	if ( this->type == "multi-dim" )
  //	{	MultiDimField *that = (MultiDimField*) this;
  //		return that->operator&& (f);                    }
	assert (0);                                                 }

MultiDimField & NumericField::join_when_first_is_block ( BlockFieldBase *first )
{//	if ( this->type == "one-dim" )
  //	{	OneDimField *that = (OneDimField*) this;
  //		return that->join_when_first_is_block ( first );  }
  //	if ( this->type == "multi-dim" )
  //	{	MultiDimField *that = (MultiDimField*) this;
  //		return that->join_when_first_is_block ( first );   }
	assert (0);                                                   }

MultiDimField & NumericField::join_when_first_is_composed ( ComposedField *first )
{//	if ( this->type == "one-dim" )
  //	{	OneDimField *that = (OneDimField*) this;
  //		return that->join_when_first_is_composed ( first );  }
  //	if ( this->type == "multi-dim" )
  //	{	MultiDimField *that = (MultiDimField*) this;
  //		return that->join_when_first_is_composed ( first );   }
	assert (0);                                                      }

MultiDimField & NumericField::join_when_first_is_multi ( MultiDimField *first )
{//	if ( this->type == "one-dim" )
  //	{	OneDimField *that = (OneDimField*) this;
  //		return that->join_when_first_is_multi ( first );  }
  //	if ( this->type == "multi-dim" )
  //	{	MultiDimField *that = (MultiDimField*) this;
  //		return that->join_when_first_is_multi ( first );   }
	assert (0);                                                   }


MultiDimField & BlockFieldBase::operator&& (NumericField &second)
{	return second.join_when_first_is_block (this);  }

MultiDimField & ComposedField::operator&& (NumericField &f)
{	return f.join_when_first_is_composed (this);  }

MultiDimField & MultiDimField::operator&& (NumericField &f)
{	return f.join_when_first_is_multi (this);  }


MultiDimField & BlockFieldBase::join_when_first_is_block (BlockFieldBase *first)
// this == second
{	ComposedField *base = new ComposedField (2);
	base->fields [0] = first;
	base->fields [1] = this;
	return * ( new MultiDimField (base) );        }

MultiDimField & BlockFieldBase::join_when_first_is_composed (ComposedField *first)
// this == second
{	size_t fvfs = first->fields.size(), n = fvfs + 1;
	ComposedField *base = new ComposedField (n);
	size_t ind = 0;
	for ( ; ind < fvfs; ind++)
	{	assert ( ind < n );
		base->fields [ind] = first->fields [ind];  }
	assert ( ind == n-1 );
	base->fields [ind] = this;
	return * ( new MultiDimField (base) );                                   }

MultiDimField & BlockFieldBase::join_when_first_is_multi (MultiDimField *first)
// this == second
{	return first->join_this_block (first, this); }


MultiDimField & ComposedField::join_when_first_is_block (BlockFieldBase *first)
// this == second
{	size_t tvfs = this->fields.size(), n = 1 + tvfs;
	ComposedField *base = new ComposedField (n);
	base->fields [0] = first;
	size_t ind = 1;
	for ( size_t i = 0; i < tvfs; i++, ind++)
	{	assert ( ind < n );
		base->fields [ind] = this->fields [i];  }
	assert ( ind == n );
	return * ( new MultiDimField (base) );                             }

MultiDimField & ComposedField::join_when_first_is_composed (ComposedField *first)
// this == second
{	size_t tvfs = this->fields.size(),
		fvfs = first->fields.size(), n = tvfs + fvfs;
	ComposedField *base = new ComposedField (n);
	size_t ind = 0;
	for ( size_t i = 0; i < fvfs; i ++, ind++)
	{	assert ( ind < n );
		base->fields [ind] = first->fields [i];  }
	for ( size_t i = 0; i < tvfs; i ++, ind++)
	{	assert ( ind < n );
		base->fields [ind] = this->fields [i];  }
	assert ( ind == n );
	return * ( new MultiDimField (base) );                                    }
	
MultiDimField & ComposedField::join_when_first_is_multi (MultiDimField *first)
// this == second
{	return first->join_this_composed (first, this); }


MultiDimField & MultiDimField::join_when_first_is_block (BlockFieldBase *first)
// this == second
{	return this->join_block_this (first, this);  }

MultiDimField & MultiDimField::join_when_first_is_composed (ComposedField *first)
// this == second
{	return this->join_composed_this (first, this);  }
	
MultiDimField & MultiDimField::join_when_first_is_multi (MultiDimField *first)
// this == second
{	return this->join_multi_this (first, this); }


///////////////////////////////////////////////////////////////
//////////////   interpolation on the Mesh   //////////////////
///////////////////////////////////////////////////////////////

// we forbid (at run time) instantiation of NumericField objects
void NumericField::interpolate ( Cell&, Cell&, double, Cell&, double )
{	assert (0);  }
void NumericField::interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double )
{	assert (0);  }
void NumericField::interpolate
	( Cell&, Cell&, double, Cell&, double, Cell&, double, Cell&, double )
{	assert (0);  }
void NumericField::interpolate
	( Cell&, Cell&, double, Cell&, double, Cell&, double, Cell&, double, Cell&, double, Cell&, double )
{	assert (0);  }


void BlockFieldBase::interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	// assert (frac_1 + frac_2 == 1.); // does not work for floating point numbers
	for ( size_t i = this->index_min; i < this->index_max_plus_one; i++ )
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i];                                 }

void BlockFieldBase::interpolate ( Cell &final, Cell &cell_1, double frac_1,
	Cell &cell_2, double frac_2, Cell &cell_3, double frac_3 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 == 1.);
	// does not work for floating point numbers
	for ( size_t i = this->index_min; i < this->index_max_plus_one; i++ )
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i];                                     }

void BlockFieldBase::interpolate ( Cell &final, Cell &cell_1, double frac_1,
	Cell &cell_2, double frac_2, Cell &cell_3, double frac_3, Cell &cell_4, double frac_4 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	assert ( (frac_4 >= 0.) && (frac_4 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 + frac_4 == 1.);
	// does not work for floating point numbers
	for ( size_t i = this->index_min; i < this->index_max_plus_one; i++ )
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i]
			+ frac_4 * cell_4.real_heap[i];                                     }

void BlockFieldBase::interpolate ( Cell &final, Cell &cell_1, double frac_1,
	Cell &cell_2, double frac_2, Cell &cell_3, double frac_3,
  Cell &cell_4, double frac_4, Cell &cell_5, double frac_5, Cell &cell_6, double frac_6 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	assert ( (frac_4 >= 0.) && (frac_4 <= 1.) );
	assert ( (frac_5 >= 0.) && (frac_5 <= 1.) );
	assert ( (frac_6 >= 0.) && (frac_6 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 + frac_4 + frac_5 + frac_6 == 1.);
	// does not work for floating point numbers
	for ( size_t i = this->index_min; i < this->index_max_plus_one; i++ )
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i]
			+ frac_4 * cell_4.real_heap[i]
			+ frac_5 * cell_5.real_heap[i]
			+ frac_6 * cell_6.real_heap[i];                                     }


void ComposedField::interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	// assert (frac_1 + frac_2 == 1.); // does not work for floating point numbers
	std::vector <BlockFieldBase*> :: iterator
		it = this->fields.begin(),
		e = this->fields.end();
	for ( ; it != e; it++ )
	{	BlockFieldBase &f = **it;
		assert (f.size() == 1);
		assert (final.dim == f.lives_on_cells_of_dim);
		assert (cell_1.dim == f.lives_on_cells_of_dim);
		assert (cell_2.dim == f.lives_on_cells_of_dim);
		size_t i = f.index_min;
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i];                  }     }

void ComposedField::interpolate
(	Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2,
	Cell &cell_3, double frac_3                                             )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 == 1.);
	// does not work for floating point numbers
	std::vector <BlockFieldBase*> :: iterator
		it = this->fields.begin(),
		e = this->fields.end();
	for ( ; it != e; it++ )
	{	BlockFieldBase &f = **it;
		assert (f.size() == 1);
		assert (final.dim == f.lives_on_cells_of_dim);
		assert (cell_1.dim == f.lives_on_cells_of_dim);
		assert (cell_2.dim == f.lives_on_cells_of_dim);
		assert (cell_3.dim == f.lives_on_cells_of_dim);
		size_t i = f.index_min;
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i];              }   }

void ComposedField::interpolate
(	Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2,
	Cell &cell_3, double frac_3, Cell &cell_4, double frac_4               )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	assert ( (frac_4 >= 0.) && (frac_4 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 + frac_4 == 1.);
	// does not work for floating point numbers
	std::vector <BlockFieldBase*> :: iterator
		it = this->fields.begin(),
		e = this->fields.end();
	for ( ; it != e; it++ )
	{	BlockFieldBase &f = **it;
		assert (f.size() == 1);
		assert (final.dim == f.lives_on_cells_of_dim);
		assert (cell_1.dim == f.lives_on_cells_of_dim);
		assert (cell_2.dim == f.lives_on_cells_of_dim);
		assert (cell_3.dim == f.lives_on_cells_of_dim);
		assert (cell_4.dim == f.lives_on_cells_of_dim);
		size_t i = f.index_min;
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i]
			+ frac_4 * cell_4.real_heap[i];                   }   }

void ComposedField::interpolate ( Cell &final, Cell &cell_1, double frac_1,
	Cell &cell_2, double frac_2, Cell &cell_3, double frac_3,
  Cell &cell_4, double frac_4, Cell &cell_5, double frac_5, Cell &cell_6, double frac_6 )
{	assert ( (frac_1 >= 0.) && (frac_1 <= 1.) );
	assert ( (frac_2 >= 0.) && (frac_2 <= 1.) );
	assert ( (frac_3 >= 0.) && (frac_3 <= 1.) );
	assert ( (frac_4 >= 0.) && (frac_4 <= 1.) );
	assert ( (frac_5 >= 0.) && (frac_5 <= 1.) );
	assert ( (frac_6 >= 0.) && (frac_6 <= 1.) );
	// assert (frac_1 + frac_2 + frac_3 + frac_4 + frac_5 + frac_6 == 1.);
	// does not work for floating point numbers
	std::vector <BlockFieldBase*> :: iterator
		it = this->fields.begin(),
		e = this->fields.end();
	for ( ; it != e; it++ )
	{	BlockFieldBase &f = **it;
		assert (f.size() == 1);
		assert (final.dim == f.lives_on_cells_of_dim);
		assert (cell_1.dim == f.lives_on_cells_of_dim);
		assert (cell_2.dim == f.lives_on_cells_of_dim);
		assert (cell_3.dim == f.lives_on_cells_of_dim);
		assert (cell_4.dim == f.lives_on_cells_of_dim);
		assert (cell_5.dim == f.lives_on_cells_of_dim);
		assert (cell_6.dim == f.lives_on_cells_of_dim);
		size_t i = f.index_min;
		final.real_heap[i] = frac_1 * cell_1.real_heap[i]
			+ frac_2 * cell_2.real_heap[i]
			+ frac_3 * cell_3.real_heap[i]
			+ frac_4 * cell_4.real_heap[i]
			+ frac_5 * cell_5.real_heap[i]
			+ frac_6 * cell_6.real_heap[i];                   }   }


void MultiDimField::interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2 )
{	MultiDimField::interpolate2 ( this, final, cell_1, frac_1, cell_2, frac_2 );  }

void MultiDimField::interpolate
(	Cell &final, Cell &cell_1, double frac_1, Cell &cell_2,
	double frac_2, Cell &cell_3, double frac_3               )
{	MultiDimField::interpolate3
		( this, final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3 );  }

void MultiDimField::interpolate
(	Cell &final, Cell &cell_1, double frac_1, Cell &cell_2,
	double frac_2, Cell &cell_3, double frac_3, Cell &cell_4, double frac_4 )
{	MultiDimField::interpolate4
		( this, final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3, cell_4, frac_4 );  }

void MultiDimField::interpolate ( Cell &final, Cell &cell_1, double frac_1,
	Cell &cell_2, double frac_2, Cell &cell_3, double frac_3,
  Cell &cell_4, double frac_4, Cell &cell_5, double frac_5, Cell &cell_6, double frac_6 )
{	MultiDimField::interpolate6 ( this, final, cell_1, frac_1, cell_2, frac_2,
		cell_3, frac_3, cell_4, frac_4, cell_5, frac_5, cell_6, frac_6 );                      }


/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//       H       H     H     H    H   H   HHHHHHH   HHH     H       HHHH       HHHH        //
//       HH     HH     H     HH   H   H   H        H   H    H       H   HH    H    H       //
//       H H   H H    H H    HH   H   H   H       H     H   H       H     H    H           //
//       H  H H  H    H H    H H  H   H   H       H     H   H       H     H     HHH        //
//       H   H   H   HHHHH   H H  H   H   HHHHH   H     H   H       H     H        H       //
//       H       H   H   H   H  H H   H   H       H     H   H       H     H         H      //
//       H       H  H     H  H   HH   H   H        H   H    H       H   HH    H   HH       //
//       H       H  H     H  H    H   H   H         HHH     HHHHHH  HHHH       HHH         //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////


// Manifolds are objects which mimick the idea of a Riemannian manifold.
// A little bit of theory : a Riemannian manifold is a differentiable manifold
// endowed with a Riemannian metric. That is, for each point of the manifold
// an inner product is given on the tangent space to that manifold at that point.
// This induces a structure of metric space (a distance) on the manifold itself.


///////////////////////////////////////////////////////////////
///////////////     factory functions     /////////////////////
///////////////////////////////////////////////////////////////


Manifold Manifold::euclid ( short int d )  // static

// trivial manifold (Euclidean space of dimension d)

{	Manifold res (d);
	res.dim = d;
	res.inner_prod_pointer = Manifold::trivial_inner_prod;
	res.distance_pointer = Manifold::trivial_distance;
	res.project_pointer = Manifold::trivial_projection;
	res.rotate_90_pointer = Manifold::trivial_rotate_90;
	return res;                                             }


Manifold & Manifold::implicit ( FunctionOnMesh::Function & psi )

// a manifold defined implicitly by the equation psi == 0
// (level set)

{	size_t n = this->coord_field->size(),
	       m = psi.components.size();
	assert ( m == 1 );
	ImplicitManifold * that = new ImplicitManifold ( n*m );
	that->outer_space = this;
	that->coord_field = this->coord_field;
	that->coord_func = this->coord_func;
	that->level_function = & psi;
	// deriv_lev_func is initialized as a vector field of dimension n
	// in the constructor of ImplicitManifold
	for ( size_t j = 0; j < m; j++ )
	for ( size_t i = 0; i < n; i++ )
		that->deriv_lev_func.components[i+j*n] =
			psi.components[j]->deriv(this->coord_func->components[i]);
	that->inner_prod_pointer = Manifold::implicit_mani_inner_prod;
	that->distance_pointer = Manifold::implicit_mani_distance;
	that->project_pointer = Manifold::implicit_mani_projection;
	that->rotate_90_pointer = Manifold::implicit_mani_rotate_90;
	Manifold & th = * ( (Manifold*) that );
	return th;                                                    }


double Manifold::trivial_inner_prod  // static

( void * that, const std::vector <double> & P, const std::vector <double> & vec1,
                                  const std::vector <double> & vec2 )
{	assert ( vec1.size() == vec2.size() );
	assert ( P.size() == vec1.size() );
	double res = 0.0;
	std::vector<double>::const_iterator it1 = vec1.begin(),
		it1_e = vec1.end(), it2 = vec2.begin();
	for ( ; it1 != it1_e; it1++, it2++ ) res += (*it1) * (*it2);
	return res;                                                    }
	

double Manifold::trivial_distance  // static
( void * that, const std::vector <double> & P, const std::vector <double> & Q )

{	assert ( P.size() == Q.size() );
	double res = 0.0;
	std::vector<double>::const_iterator it1 = P.begin(),
		it1_e = P.end(), it2 = Q.begin();
	for ( ; it1 != it1_e; it1++, it2++ )
	{	double dif = *it2 - *it1;
		res += dif * dif;           }
	return sqrt(res);                                     }


std::vector <double> & Manifold::trivial_projection  // static
( void * that, const std::vector <double> & x )

{	std::vector <double> *y = new std::vector <double> (x);
	return *y;                                                  }


std::vector <double> & Manifold::trivial_rotate_90  // static
( void * that, const std::vector <double> & point, const std::vector <double> & vec, double desired_norm )

{	assert ( point.size() == 2 );
	assert ( vec.size() == 2 );
	std::vector <double> * res = new std::vector<double> (2);
	(*res)[0] = vec[1];  (*res)[1] = -vec[0];
	double norm = 0.0;
	std::vector<double>::iterator it = res->begin(), it_e = res->end();
	for ( ; it != it_e; it++ ) norm += (*it) * (*it);
	norm = sqrt ( norm );
	double frac = desired_norm / norm;
	it = res->begin(), it_e = res->end();
	for ( ; it != it_e; it++ ) *it *= frac;
	return * res;                                                       }


double Manifold::implicit_mani_inner_prod  // static
( void * that, const std::vector <double> & P, const std::vector <double> & vec1,
                                  const std::vector <double> & vec2 )
// copied from trivial
{	assert ( vec1.size() == vec2.size() );
	assert ( P.size() == vec1.size() );
	double res = 0.0;
	std::vector<double>::const_iterator it1 = vec1.begin(),
		it1_e = vec1.end(), it2 = vec2.begin();
	for ( ; it1 != it1_e; it1++, it2++ ) res += (*it1) * (*it2);
	return res;                                                    }
	

double Manifold::implicit_mani_distance  // static
( void * that, const std::vector <double> & P, const std::vector <double> & Q )
// copied from trivial
{	assert ( P.size() == Q.size() );
	double res = 0.0;
	std::vector<double>::const_iterator it1 = P.begin(),
		it1_e = P.end(), it2 = Q.begin();
	for ( ; it1 != it1_e; it1++, it2++ )
	{	double dif = *it2 - *it1;
		res += dif * dif;           }
	return sqrt(res);                                     }


std::vector <double> & Manifold::implicit_mani_projection  // static
( void * that, const std::vector <double> & x )

{	ImplicitManifold & manif = * ( (ImplicitManifold*) that );
	std::vector <double> & y = * ( new std::vector <double> (x) );
	// we project using Newton's method for undetermined systems
	FunctionOnMesh::Function & psi = * ( manif.level_function );
	assert ( psi.components.size() == 1 );
	size_t n = manif.coord_field->size();
	assert ( n == x.size() );
	FunctionOnMesh::Function & coord = * ( manif.coord_func );
	Cell & P = Cell::point();
	for ( short int k = 0; k < ImplicitManifold::steps_for_Newton; k++ )
	{	for ( size_t i = 0; i < n; i++ ) coord[i] == y[i];  // at P
		double psi_at_P = psi(P);
		std::vector <double> gradpsi_at_P(n);
		for ( size_t i = 0; i < n; i++ )
			gradpsi_at_P[i] = manif.deriv_lev_func[i](P);
		double norm2 = 0.;
		for ( size_t i = 0; i < n; i++ )
			norm2 += gradpsi_at_P[i]*gradpsi_at_P[i];
		double coef = psi(P) / norm2;
		for ( size_t i = 0; i < n; i++ )
			y[i] -= coef * gradpsi_at_P[i];                        }
	return y;                                                                   }


std::vector <double> & Manifold::implicit_mani_rotate_90  // static
( void * data, const std::vector <double> & point, const std::vector <double> & vec, double desired_norm )
// copied from trivial
{	assert ( point.size() == 2 );
	assert ( vec.size() == 2 );
	std::vector <double> * res = new std::vector<double> (2);
	(*res)[0] = vec[1];  (*res)[1] = -vec[0];
	double norm = 0.0;
	std::vector<double>::iterator it = res->begin(), it_e = res->end();
	for ( ; it != it_e; it++ ) norm += (*it) * (*it);
	norm = sqrt ( norm );
	double frac = desired_norm / norm;
	it = res->begin(), it_e = res->end();
	for ( ; it != it_e; it++ ) *it *= frac;
	return * res;                                                       }


FunctionOnMesh::Function & Manifold::coordinate_system ( std::string s )

// builds a FunctionOnMesh, with a MultiDimField associated
// 's' may be "Lagrange degree one" for instance

{	assert ( s == "Lagrange degree one" );
	assert ( this->dim > 1 );
	MultiDimField & fi = NumericField::multi_dim ("size", this->dim, "lives on", "points");
	this->coord_field = &fi;
	FunctionOnMesh::Function & fu = FunctionOnMesh::from_field ( fi, "Lagrange degree one" );
	this->coord_func = &fu;
	Cell::init_cell[0]->push_back ( & Mesh::prescribe_on );
	Cell::data_for_init[0]->push_back ( NULL );
	// what if we call this more than once ? damn ...
	return fu;                                                                                }

