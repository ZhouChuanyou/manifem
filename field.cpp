
// maniFEM field.cpp 2019.10.30

#include "field.h"

using namespace maniFEM;

size_t Field::Scalar::nb_of_components ( )  // virtual from Field::Core
{	return 1;  }


size_t Field::Block::nb_of_components ( )  // virtual from Field::Corex
{	return Mesh::diff ( this->max_index_p1, this->min_index );  }
// Mesh::diff  provides a safe way to substract two size_t numbers

Field::Scalar * Field::Scalar::component ( size_t i )  // virtual from Field::Core
{	assert ( i == 0 );
	return this;        }

Field::Scalar * Field::Block::component ( size_t i )  // virtual from Field::Core
// we build a new Field::Scalar
// it would be nice to check whether a field for index i exists already ...
{	assert ( i < this->nb_of_components() );
	return new Field::Scalar ( tag::lives_on_positive_cells, tag::of_dim,
	  this->lives_on_cells_of_dim, tag::has_index_in_heap, this->min_index + i );  }





