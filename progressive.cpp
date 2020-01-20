
// manifem/progressive.cpp 2020.01.13

#include <stack>
#include "math.h"
#include <random>

#include "maniFEM.h"
#include "progressive.h"

using namespace maniFEM;


std::map < Cell::Core *, Mesh::Officer::RankZero * > Mesh::Officer::officer;  // use Cell::hook

double Mesh::Cloud::ratio_between_distances; // must be > 2
std::vector <double> Mesh::Cloud::distance;
std::vector <double> Mesh::Cloud::range;
Mesh::Officer * Mesh::Cloud::highest_officer;

std::map < Cell::Core *, std::vector < double > > Mesh::Progressive::normals;  // use Cell::hook

Cell Mesh::Progressive::temporary_vertex ( tag::non_existent );

//-------------------------------------------------------------------------------------------------

void Mesh::Cloud::add_point ( Cell p )

{	assert ( p.dim() == 0 );
	assert ( Mesh::Cloud::distance.size() );
	assert ( Mesh::Cloud::distance.size() == Cloud::range.size() );
	Mesh::Officer * of = new Mesh::Officer::RankZero ( p );

	if ( Mesh::Cloud::highest_officer == nullptr )  // cloud previously empty
	{	assert ( of->superior == nullptr );
		Mesh::Cloud::highest_officer = of;
		return;                              }

	if ( Mesh::Cloud::highest_officer->get_rank() == 0 )
	{	Mesh::Cloud::build_vectors ( 1 );
		Mesh::Cloud::highest_officer = Mesh::Cloud::highest_officer->promote ( 1 );  }
	// perhaps we should define 'adopt' as a method of class Mesh::Officer
	// this way, we could eliminate the above check
	
	assert ( Mesh::Cloud::highest_officer->get_rank() > 0 );
	Mesh::Officer::PositiveRank * hof =
		static_cast < Mesh::Officer::PositiveRank * > ( Mesh::Cloud::highest_officer );
	hof->adopt ( of );
	
}  // end of Cloud::add_point

//-------------------------------------------------------------------------------------------------

void Mesh::Officer::PositiveRank::adopt ( Mesh::Officer * of, Mesh::Officer * to_avoid )

// 'to_avoid' defaults to nullptr, so may be called with one argument only

// perhaps we should define 'adopt' as a method of class Mesh::Officer
// this way, we could eliminate the check in 'add_point'
// also, define 'rank' as a data member in Mesh::Officer

{	size_t rank_of_of = of->get_rank();
	size_t rank_of_this = this->get_rank();
	assert ( rank_of_this >= rank_of_of );
	double d = Manifold::working.distance ( of->point, this->point );

	if ( rank_of_this > rank_of_of )
		if ( d <= Mesh::Cloud::range[rank_of_this] )
		{	assert ( rank_of_this );
			Mesh::Officer::PositiveRank * hof =
				static_cast < Mesh::Officer::PositiveRank * > ( this );
			bool success = of->try_to_cling_to ( hof, to_avoid );
			if ( success ) return;                                      }

	// 'of' is out of the range of 'this'
	// however, 'this' may have a superior
	if ( this->superior )
	{	Mesh::Officer::PositiveRank * hof = this->superior;
		hof->adopt ( of, this );  // avoid 'this' when trying to cling to subalterns
		return;                                               }

	// 'of' is out of the range of the Cloud
	assert ( this == Mesh::Cloud::highest_officer );  // and we want to promote 'this'
	size_t size_of_dist = Mesh::Cloud::distance.size();
	assert ( size_of_dist > rank_of_this );
	assert ( Mesh::Cloud::range.size() == size_of_dist );
	double dd = Mesh::Cloud::distance[rank_of_this];
	while ( true )  // first, we use 'distance' and 'range' at their existing size
	{	rank_of_this++;
		if ( rank_of_this >= size_of_dist ) break;
		dd = Mesh::Cloud::distance[rank_of_this];
		if ( d <= dd ) break;                       }
	if ( d > dd )  // we need to promote more, while increasing the size of 'distance' and 'range'
		while ( true )
		{	assert ( Mesh::Cloud::distance.size() == rank_of_this );
			assert ( Mesh::Cloud::range.size() == rank_of_this );
			dd *= Mesh::Cloud::ratio_between_distances;
			Mesh::Cloud::distance.push_back ( dd );
			Mesh::Cloud::range.push_back ( Mesh::Cloud::range.back() + dd );
			if ( d <= dd ) break;
			rank_of_this++;                                                   }

	assert ( Mesh::Cloud::distance.size() >= rank_of_this + 1 );
	assert ( Mesh::Cloud::distance.size() == Mesh::Cloud::range.size() );
	Mesh::Officer * new_this = this->promote ( rank_of_this );
	// 'promote' builds the whole chain of intermediate officers
	assert ( new_this->get_rank() == rank_of_this );
	assert ( new_this->superior == nullptr );
	Mesh::Cloud::highest_officer = new_this;
	assert ( new_this->get_rank() );
	assert ( d <= Mesh::Cloud::distance[rank_of_this] );
	Mesh::Officer::PositiveRank * hof =
		static_cast < Mesh::Officer::PositiveRank * > ( new_this );
	hof->direct_adopt ( of, rank_of_this );
	
}  // end of Officer::PositiveRank::adopt

//-------------------------------------------------------------------------------------------------

bool Mesh::Officer::try_to_cling_to
( Mesh::Officer::PositiveRank * insert_place, Mesh::Officer * to_avoid )

// 'to_avoid' may be nullptr, meaning there is nothing to avoid
	
{	size_t rank_of_this = this->get_rank();
	size_t rank_of_ip = insert_place->get_rank();
	assert ( Mesh::Cloud::distance.size() > rank_of_ip );
	assert ( rank_of_ip > rank_of_this );

	assert ( Manifold::working.distance ( this->point, insert_place->point )
	  < Mesh::Cloud::range[rank_of_ip] );

	// try to cling to subalterns
	// but not to all of them, only to those who are not too far
	// also, there might be one of them not worth analysing (argument to_avoid)
	size_t rank_of_sub = rank_of_ip - 1;
	if ( rank_of_sub > rank_of_this )
	{	std::set < Mesh::Officer * > & ipsub = insert_place->subalterns;
		std::set<Mesh::Officer*>::iterator it = ipsub.begin();
		for ( ; it != ipsub.end(); it++ )
		{	Mesh::Officer * sub = *it;
			assert ( sub );
			if ( sub == to_avoid ) continue;
			assert ( rank_of_sub == sub->get_rank() );
			if ( Manifold::working.distance ( this->point, sub->point )
			  > Mesh::Cloud::range[rank_of_sub] ) continue;
			assert ( rank_of_sub );
			Mesh::Officer::PositiveRank * sb =
				static_cast < Mesh::Officer::PositiveRank * > ( sub );
			bool success = this->try_to_cling_to ( sb, nullptr );  // nothing to avoid
			if ( success ) return true;                                                        }  }

	// couldn't cling to any subaltern; try 'insert_place' itself
	if ( Manifold::working.distance ( this->point, insert_place->point )
			 > Mesh::Cloud::distance[rank_of_ip] )
		return false;  // no success
	insert_place->direct_adopt ( this, rank_of_ip );
	return true;  // success

}  // end of Officer::try_to_cling_to

//-------------------------------------------------------------------------------------------------

void Mesh::Cloud::remove_point ( Cell p )

{	assert ( p.dim() == 0 );
	Mesh::Officer * of_0 = Mesh::Officer::officer[p.core];
// check ! but we will use Cell::hook instead ...
	Mesh::Officer::officer.erase ( p.core );
	assert ( of_0 );
	assert ( of_0->point.core == p.core );
	
	// look for higher rank officers exactly above this point
	Mesh::Officer::PositiveRank * of_1 = of_0->superior;
	if ( of_1 == nullptr )  // of_0 is the last and only point in the cloud
	{	delete ( of_0 );
		assert ( Mesh::Cloud::highest_officer == of_0 );
		Mesh::Cloud::highest_officer = nullptr;  // cloud becomes empty
		return;                                    }

	// begin to remove officers, from bottom to top
	while ( of_1 )
	{	size_t nn = of_1->subalterns.erase ( of_0 );
		assert ( nn == 1 );
		if ( of_0->get_rank() )  // 'of_1' must adopt subalterns of 'of_0'
		{	Mesh::Officer::PositiveRank * of_00 =
				static_cast < Mesh::Officer::PositiveRank * > ( of_0 );
			std::set<Mesh::Officer*>::iterator it = of_00->subalterns.begin();
			for ( ; it != of_00->subalterns.end(); it++ )
			{	assert ( *it );  of_1->adopt ( *it );  }                     }
		delete of_0;
		of_0 = of_1;
		if ( of_0->point.core != p.core ) break;
		of_1 = of_1->superior;                                                }

	if ( of_1 == nullptr )  // normal exit from the above loop
	{	assert ( of_0 == Mesh::Cloud::highest_officer );
		assert ( of_0->point.core == p.core );
		// 'of_0' will be deleted, one of its subalterns will adopt the others
		assert ( of_0->get_rank() );
		of_1 = static_cast < Mesh::Officer::PositiveRank * > ( of_0 );
		std::set<Mesh::Officer*>::iterator it = of_1->subalterns.begin();
		Mesh::Cloud::highest_officer = *it;
		Mesh::Cloud::highest_officer->superior = nullptr;
		for ( it++; it != of_1->subalterns.end(); it++ )
		{	// Cloud::highest_officer may change after each adoption
			// it will surely change after the first adoption, must get promoted
			Mesh::Officer::PositiveRank * hof =
				static_cast < Mesh::Officer::PositiveRank * > ( Mesh::Cloud::highest_officer );
		  hof->adopt ( *it );                                                                }
		of_0 = Mesh::Cloud::highest_officer;  // just for the 'if' below to be affirmative
		delete of_1;                                                                            }

	// collapse highest officer(s) if they have only one subaltern
	if ( of_0 == Mesh::Cloud::highest_officer )
	{	while ( of_0->get_rank() )
		{	of_1 = static_cast < Mesh::Officer::PositiveRank * > ( of_0 );
			if ( of_1->subalterns.size() != 1 ) break;
			of_0 = * ( of_1->subalterns.begin() );
			assert ( of_0 );
			delete of_1;
			if ( of_0->get_rank() == 0 ) break;                                }
		of_0->superior = nullptr;
		Mesh::Cloud::highest_officer = of_0;                                   }

} // end of Mesh::Cloud::remove_point

//-------------------------------------------------------------------------------------------------

void Mesh::Officer::PositiveRank::find_close_vertices  // virtual from Mesh::Officer	
( Cell cll, double dist_high, std::set < Cell::Core* > & res )

// among subalterns of 'this', returns all vertices closer to cll than dist_high
// surely returns cll itself (if it belongs to the cloud)
// and most likely also its two immediate neighbours

// 'cll' is just a cell, does not need to belong to the cloud
	
{	double d = Manifold::working.distance ( cll, this->point );
	// try to avoid use of square root !
	if ( d > Mesh::Cloud::range [ this->get_rank() ] + dist_high ) return;
	std::set<Mesh::Officer*>::iterator it = this->subalterns.begin();
	for ( ; it != this->subalterns.end(); it++ )
		(*it)->find_close_vertices ( cll, dist_high, res );               }
	
//-------------------------------------------------------------------------------------------------

void Mesh::Officer::RankZero::find_close_vertices  // virtual from Mesh::Officer	
( Cell cll, double dist_high, std::set < Cell::Core* > & res )

// just insert itself in the set of results, if close enough to cll
	
{	double d = Manifold::working.distance ( cll, this->point );
	if ( d > dist_high ) return;
	res.insert ( this->point.core );                                  }
	
//-------------------------------------------------------------------------------------------------

void Mesh::Cloud::print_everything ()
{
	std::cout << "ratio between distances : " << Mesh::Cloud::ratio_between_distances << std::endl << std::flush;
	std::cout << "rank distances :";
	std::vector<double>::iterator it = Mesh::Cloud::distance.begin(),
		it_e = Mesh::Cloud::distance.end();
	for ( ; it != it_e; it++ )  std::cout << " " << *it;
	std::cout << std::endl;
	std::cout << "ranges :";
	it = Mesh::Cloud::range.begin();
	it_e = Mesh::Cloud::range.end();
	for ( ; it != it_e; it++ )  std::cout << " " << *it;
	std::cout << std::endl;
	if ( Mesh::Cloud::highest_officer )
	{	if ( Mesh::Cloud::highest_officer->get_rank() > 0 )
		{	std::cout << "officer hierarchy : " << std::endl;
			Mesh::Officer::PositiveRank * hro =
				static_cast < Mesh::Officer::PositiveRank * > ( Mesh::Cloud::highest_officer );
			hro->print_subalterns ( 0 );                                               }
		else std::cout << "only one rank zero officer, " << Mesh::Cloud::highest_officer->point.core->name << std::endl;  }
	else std::cout << "empty cloud" << std::endl;
}

//-------------------------------------------------------------------------------------------------

void Mesh::Officer::PositiveRank::print_subalterns ( size_t indent_level )
{
	size_t local_rank = this->get_rank();
	for ( size_t i = 0; i < indent_level; i++ )  std::cout << "  ";
	std::cout << this->point.core->name << "_r" << rank << " has " << this->subalterns.size()
            << " subaltern";
	if ( this->subalterns.size() > 1 ) std::cout << "s";
	std::cout << std::endl << std::flush;	
	if ( local_rank == 1 )
	{	for ( size_t i=0; i <= indent_level; i++ )  std::cout << "  ";
		std::set<Mesh::Officer*>::iterator it = this->subalterns.begin(),
			it_e = this->subalterns.end();
		for ( ; it != it_e; it++ )
		{	Cell p = (*it)->point;
			std::cout << p.core->name << " " ;  }
		std::cout << std::endl << std::flush;                 }
	else // local_rank > 1
	{	indent_level++;
		std::set<Mesh::Officer*>::iterator it = this->subalterns.begin(),
			it_e = this->subalterns.end();
		for ( ; it != it_e; it++ )
		{	Mesh::Officer::PositiveRank * of =
				static_cast < Mesh::Officer::PositiveRank * > ( *it );
			of->print_subalterns ( indent_level );                   }   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

inline void improve_normal  // static
( const Cell & A, std::vector < double > & nor )

// project 'nor' onto working manifold

{	// we normalize 'nor'
	double n2 = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		n2 += nor[i]*nor[i];  // use Riemannian metric !
	double norm = approx_sqrt ( n2,
	  tag::around, Mesh::Progressive::desired_length_sq, Mesh::Progressive::desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = Mesh::Progressive::desired_length / norm;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  nor[i] *= norm;
	// we project on Manifold::working
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( Mesh::Progressive::temporary_vertex ) = x(A) + nor[i];  }
	Manifold::working.project ( Mesh::Progressive::temporary_vertex );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		nor[i] = x ( Mesh::Progressive::temporary_vertex ) - x(A);  }
}

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::improve_normal  // static
( const Cell & A, const Cell & B, std::vector < double > & AB_coord, std::vector < double > & nor )

// project 'nor' onto working manifold and normalize it relatively to segment AB

{	// we normalize 'nor'
	double n2 = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		n2 += nor[i]*nor[i];  // use Riemannian metric !
	double norm = approx_sqrt ( n2,
	  tag::around, Mesh::Progressive::desired_length_sq, Mesh::Progressive::desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = Mesh::Progressive::desired_length / norm;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  nor[i] *= norm;
	// we project on Manifold::working
	std::vector < double > mid_seg ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		mid_seg[i] = ( x(A) + x(B) ) / 2.;
		x ( Mesh::Progressive::temporary_vertex ) = mid_seg[i] + nor[i];  }
	Manifold::working.project ( Mesh::Progressive::temporary_vertex );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		nor[i] = Manifold::working.coordinates()[i] ( Mesh::Progressive::temporary_vertex ) - mid_seg[i];
	// we ensure again that 'nor' is orthogonal to AB_coord
	double prod = 0.;  n2 = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  // use Riemannian metric !
	{	prod += nor[i] * AB_coord[i];  n2 += AB_coord[i] * AB_coord[i];    }
	prod /= n2;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ ) nor[i] -= prod * AB_coord[i];
	// we normalize 'nor'
	n2 = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  n2 += nor[i]*nor[i];  // use Riemannian metric !
	norm = approx_sqrt ( n2,
	  tag::around, Mesh::Progressive::desired_length_sq, Mesh::Progressive::desired_length );
	norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
	norm = Mesh::Progressive::desired_length / norm;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  nor[i] *= norm;
}

//-------------------------------------------------------------------------------------------------

inline bool opposite_signs ( const double a, const double b )
{	if ( a < 0. )  return b >= 0.;
	if ( a == 0. )  return true;
	return b <= 0.;                  }

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver_c1 ( const double & desired_length )

// search for a starting point in a manifold of co-dimension one
// e.g. a curve in the plane or a surface in 3D

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	// Function m_impl->level_function;

	const size_t nc = Manifold::working.coordinates().nb_of_components();
	const double len_sq = desired_length * desired_length;
	
	Cell tmp_ver ( tag::vertex );
	Cell tmp_ver_1 ( tag::vertex );
	Cell tmp_ver_2 ( tag::vertex );

	std::default_random_engine random_generator;
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_length;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		for ( size_t j = 0; j < nc; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_ver_1 ) = distr(random_generator);
				x ( tmp_ver_2 ) = distr(random_generator);        }
			double v1 = m_impl->level_function ( tmp_ver_1 ),
			       v2 = m_impl->level_function ( tmp_ver_2 );
			if ( opposite_signs ( v1, v2 ) )
				// refine by applying bissection algorithm
				while ( true )
				{	if ( Manifold::working.dist_sq ( tmp_ver_1, tmp_ver_2 ) < len_sq )
					{	tmp_ver.dispose();  tmp_ver_2.dispose();
						Manifold::working.project ( tmp_ver_1 );
						return tmp_ver_1;                          }
					m_impl->surrounding_space.interpolate ( tmp_ver, 0.5, tmp_ver_1, 0.5, tmp_ver_2 );
					double v = m_impl->level_function ( tmp_ver );
					if ( opposite_signs ( v, v2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_1 ) = x ( tmp_ver );                 }
						v1 = v;                                               }
					else if ( opposite_signs ( v, v1 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_ver_2 ) = x ( tmp_ver );                 }
						v2 = v;                                               }
					else assert( false );                                                }  }
		size_of_cube *= 2;                                                                                 }
	tmp_ver_1.dispose();  tmp_ver_2.dispose();
}

//-------------------------------------------------------------------------------------------------

inline double ext_prod_R2 ( const double vx, const double vy, const double wx, const double wy )
{	return vx*wy - wx*vy;  }


inline bool origin_outside ( const double & Ax, const double & Ay,
                             const double & Bx, const double & By,
                             const double & Cx, const double & Cy )

{
//	std::cout << "(" << Ax << "," << Ay << ") (" << Bx << "," << By << ") (" << Cx << "," << Cy << ")" << std::endl;
//	std::cout << "AB  " << ext_prod_R2 ( Bx-Ax, By-Ay, Ax, Ay ) << " " << ext_prod_R2 ( Bx-Ax, By-Ay, Cx-Ax, Cy-Ay ) << std::endl;
//	std::cout << "BC  " << ext_prod_R2 ( Cx-Bx, Cy-By, Bx, By ) << " " << ext_prod_R2 ( Cx-Bx, Cy-By, Ax-Bx, Ay-By ) << std::endl;
//	std::cout << "CA  " << ext_prod_R2 ( Ax-Cx, Ay-Cy, Cx, Cy ) << " " << ext_prod_R2 ( Ax-Cx, Ay-Cy, Bx-Cx, By-Cy ) << std::endl;
	return
	   opposite_signs ( - ext_prod_R2 ( Bx-Ax, By-Ay, Ax, Ay ), ext_prod_R2 ( Bx-Ax, By-Ay, Cx-Ax, Cy-Ay ) )
	or opposite_signs ( - ext_prod_R2 ( Cx-Bx, Cy-By, Bx, By ), ext_prod_R2 ( Cx-Bx, Cy-By, Ax-Bx, Ay-By ) )
	or opposite_signs ( - ext_prod_R2 ( Ax-Cx, Ay-Cy, Cx, Cy ), ext_prod_R2 ( Ax-Cx, Ay-Cy, Bx-Cx, By-Cy ) ); }

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver_c2 ( const double & desired_length )

// search for a starting point in a manifold of co-dimension two (a curve in 3D)

{	Manifold::Implicit::TwoEquations * m_impl =
		dynamic_cast<Manifold::Implicit::TwoEquations*> ( Manifold::working.core );
	assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );

	const size_t nc = Manifold::working.coordinates().nb_of_components();
	const double len_sq = desired_length * desired_length;
	
	Cell tmp_A ( tag::vertex );
	Cell tmp_B ( tag::vertex );
	Cell tmp_C ( tag::vertex );
	Cell tmp_AB ( tag::vertex );
	Cell tmp_BC ( tag::vertex );
	Cell tmp_CA ( tag::vertex );

	std::default_random_engine random_generator;
	size_t counter = 0;
	size_t size_of_cube = 5;
	while ( true )
	{	double s = size_of_cube * desired_length;
		std::uniform_real_distribution<double> distr ( -s, s );
		size_t nb = 1;
		restart :
		for ( size_t j = 0; j < nc; j++ ) nb *= size_of_cube;
		for ( size_t j = 0; j < nb; j++ )
		{	for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x ( tmp_A ) = distr(random_generator);
				x ( tmp_B ) = distr(random_generator);
				x ( tmp_C ) = distr(random_generator);             }
			counter++;
			double vA1 = m_impl->level_function_1 ( tmp_A ),
			       vA2 = m_impl->level_function_2 ( tmp_A );
			double vB1 = m_impl->level_function_1 ( tmp_B ),
			       vB2 = m_impl->level_function_2 ( tmp_B );
			double vC1 = m_impl->level_function_1 ( tmp_C ),
			       vC2 = m_impl->level_function_2 ( tmp_C );
			if ( not origin_outside ( vA1, vA2, vB1, vB2, vC1, vC2 ) )
				// refine by repeatedly cutting the triangle
				while ( true )
				{	if ( ( Manifold::working.dist_sq ( tmp_A, tmp_B ) < len_sq ) and
					     ( Manifold::working.dist_sq ( tmp_B, tmp_C ) < len_sq ) and
					     ( Manifold::working.dist_sq ( tmp_C, tmp_A ) < len_sq )      )
					{	tmp_B.dispose();  tmp_C.dispose();
						tmp_AB.dispose();  tmp_BC.dispose();  tmp_CA.dispose();
						Manifold::working.project ( tmp_A );
						return tmp_A;                                            }
					m_impl->surrounding_space.interpolate ( tmp_AB, 0.5, tmp_A, 0.5, tmp_B );
					double vAB1 = m_impl->level_function_1 ( tmp_AB ),
					       vAB2 = m_impl->level_function_2 ( tmp_AB );
				  m_impl->surrounding_space.interpolate ( tmp_BC, 0.5, tmp_B, 0.5, tmp_C );
					double vBC1 = m_impl->level_function_1 ( tmp_BC ),
					       vBC2 = m_impl->level_function_2 ( tmp_BC );
				  m_impl->surrounding_space.interpolate ( tmp_CA, 0.5, tmp_C, 0.5, tmp_A );
					double vCA1 = m_impl->level_function_1 ( tmp_CA ),
					       vCA2 = m_impl->level_function_2 ( tmp_CA );
					if ( not origin_outside ( vA1, vA2, vAB1, vAB2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_B ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_CA );  }
						vB1 = vAB1;  vB2 = vAB2;  vC1 = vCA1;  vC2 = vCA2;            }
					else if ( not origin_outside ( vB1, vB2, vAB1, vAB2, vBC1, vBC2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_AB );  x ( tmp_C ) = x ( tmp_BC );  }
						vA1 = vAB1;  vA2 = vAB2;  vC1 = vBC1;  vC2 = vBC2;            }
					else if ( not origin_outside ( vC1, vC2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_CA );  x ( tmp_B ) = x ( tmp_BC );  }
						vA1 = vCA1;  vA2 = vCA2;  vB1 = vBC1;  vB2 = vBC2;            }
					else if ( not origin_outside ( vAB1, vAB2, vBC1, vBC2, vCA1, vCA2 ) )
					{	for ( size_t i = 0; i < nc; i++ )
						{	Function x = Manifold::working.coordinates()[i];
							x ( tmp_A ) = x ( tmp_BC );  x ( tmp_B ) = x ( tmp_CA );
							x ( tmp_C ) = x ( tmp_AB );                               }
						vA1 = vBC1;  vA2 = vBC2;  vB1 = vCA1;  vB2 = vCA2;
						vC1 = vAB1;  vC2 = vAB2;                                                 }
					else  // nasty nonlinear level functions ...
						goto restart;                                                               }  }
		size_of_cube *= 2;                                                                                 }
}

//-------------------------------------------------------------------------------------------------

inline Cell search_start_ver ( const double & desired_length )

// search for a starting point
// current working manifold may have co-dimension one
// e.g. a curve in the plane or a surface in 3D
// of co-dimension two (a curve in 3D)

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	if ( m_impl ) return search_start_ver_c1 ( desired_length );  // co-dimension one
	return search_start_ver_c2 ( desired_length );  // co-dimension two
}

//-------------------------------------------------------------------------------------------------

inline void redistribute_vertices ( const Mesh & msh,
  const Cell & start, const Cell & stop, double desired_length, double last_length, size_t n )

// chain of n segments, desired length d, last length d'
// move p1 with   (d'-d) / n
//    alpha x0 + beta x2 == x1 + (d'-d) / n   alpha+beta == 1
//    x2 == x0 + 2d   x1 == x0 + d
//    2 d beta == d + (d'-d)/ n 
// move p2 with 2*(d'-d) / n
//    alpha x1 + beta x3 == x2 + 2(d'-d)/n   alpha+beta == 1
//    x1 == x0 + d + (d'-d)/n
//    x3 == x0 + 3d
//    alpha (d'-d)/n + 2 d beta == d + 2(d'-d)/n
//    ( 2d - (d'-d)/n ) beta == d + (d'-d)/n
// move p2 with 3*(d'-d) / n
//    ( 2d - 2(d'-d)/n ) beta == d + (d'-d)/n

{	Cell A = msh.cell_behind(stop).base().reverse();
	for ( size_t i = 1; i < n; i++ )
	{	if ( A == start )  {  n = i;  break;  }
		A = msh.cell_behind(A).base().reverse();  }
	assert ( n > 1 );
	Cell B = msh.cell_in_front_of(A).tip();
	Cell C = msh.cell_in_front_of(B).tip();
	double epsilon = ( last_length - desired_length ) / n;
	double v1 = desired_length + epsilon,  v2 = 2*desired_length;
	while ( C != stop )
	{	double beta = v1/v2;
		Manifold::working.interpolate ( B, 1.-beta, A, beta, C );
		v2 -= epsilon;  assert ( v2 > 0. );
		A = B;  B = C;
		C = msh.cell_in_front_of(B).tip();                            }
	Manifold::working.interpolate ( B, 0.5, A, 0.5, C );               }

//-------------------------------------------------------------------------------------------------

inline double get_z_baric ( const Cell & tri )

{	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	CellIterator it = tri.boundary().iter_over ( tag::vertices );
	Function z = Manifold::working.coordinates()[2];
	double zz = 0.;
	size_t counter = 0;
	for ( it.reset(); it.in_range(); it++, counter++ )  zz += z(*it);
	assert ( counter == 3 );
	return  zz/3.;                                                               }

//-------------------------------------------------------------------------------------------------

bool tri_correctly_oriented ( const Cell & tri )

{	assert ( tri.dim() == 2 );

	CellIterator it = tri.boundary().iter_over ( tag::segments );
	it.reset();  assert ( it.in_range() );
	Cell AB = *it;
	it++;  assert ( it.in_range() );
	Cell BC = *it;
	it++;  assert ( it.in_range() );
	// Cell CA = *it;
	it++;  assert ( not it.in_range() );

	Cell A = AB.base().reverse();
	Cell B = AB.tip();
	assert ( B == BC.base().reverse() );
	Cell C = BC.tip();

	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	double  xAB = x(B) - x(A),  yAB = y(B) - y(A),
	        xBC = x(C) - x(B),  yBC = y(C) - y(B);
	return  xAB * yBC > yAB * xBC;;                                  }

//-------------------------------------------------------------------------------------------------

bool correctly_oriented ( const Mesh msh )

// tells whether 'msh's orientation is consistent with the orientation of the
// surrounding Euclidian space

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	assert ( m_euclid );
	const size_t nc = Manifold::working.coordinates().nb_of_components();

	// for surfaces, we should search the vertex with zmax and check the orientation
	// of all surrounding triangles (we need an iterator over cells above)
	if ( msh.dim() != 1 )
	{	assert ( msh.dim() == 2 );
		assert ( nc == 3 );
		CellIterator it = msh.iter_over ( tag::cells_of_dim, 2 );
		it.reset();  assert ( it.in_range() );
		Cell trimax = *it;
		double zmax = get_z_baric ( trimax );
		for ( it++; it.in_range(); it++ )
		{	double zz = get_z_baric ( *it );
			if ( zz > zmax )
			{	zmax = zz;  trimax = *it;  }   }
		return tri_correctly_oriented ( trimax );                   }

	assert ( msh.dim() == 1 );
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];

	CellIterator it = msh.iter_over ( tag::vertices );
	it.reset();  assert ( it.in_range() );
	Cell ver = *it;
	double ymax = y(ver);
	for ( it++; it.in_range(); it++ )
	{	Cell other_ver = *it;
		double other_y = y(other_ver);
		if ( other_y > ymax )
		{	ymax = other_y;  ver = other_ver;  }  }
	Cell prev_seg = msh.cell_behind ( ver );
	Cell next_seg = msh.cell_in_front_of ( ver );
	assert ( prev_seg.tip() == ver );
	assert ( next_seg.base().reverse() == ver );
	Cell A = prev_seg.base().reverse();
	Cell C = next_seg.tip();
	bool prev_orient = ( x(ver) < x(A) );
	bool next_orient = ( x(C) < x(ver) );
	if ( prev_orient == next_orient )  return prev_orient;
	double  xAB = x(ver) - x(A),  yAB = y(ver) - y(A),
	        xBC = x(C) - x(ver),  yBC = y(C) - y(ver);
	return  xAB * yBC > yAB * xBC;;
	
}

//-------------------------------------------------------------------------------------------------

bool correctly_oriented_complicated ( const Mesh msh )

// tells whether 'msh's orientation is consistent with the orientation of the
// surrounding Euclidian space

{	Manifold::Implicit::OneEquation * m_impl =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	assert ( m_impl );
	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	assert ( m_euclid );
	const size_t nc = Manifold::working.coordinates().nb_of_components();

	if ( nc != 2 )
	{	std::cout << "for the moment I can only check the orientation of "
	            << "(closed) curves in the plane - sorry" << std::endl;
		exit ( 1 );                                                         }

	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	CellIterator it = msh.iter_over ( tag::segments );
	it.reset();  assert ( it.in_range() );
	Cell seg = *it;
	assert ( seg.dim() == 1 );
	bool dx_pos = ( x (seg.tip()) - x (seg.base().reverse()) ) > 0.;
	bool dy_pos = ( y (seg.tip()) - y (seg.base().reverse()) ) > 0.;
	int counter = 0;
	for ( it++; it.in_range(); it++ )
	{	seg = *it;
		bool dx_pos_now = ( x(seg.tip()) - x (seg.base().reverse()) ) > 0.;
		bool dy_pos_now = ( y(seg.tip()) - y (seg.base().reverse()) ) > 0.;
		if ( dx_pos_now != dx_pos )  //  dx has changed sign
		{	if ( dy_pos_now != dy_pos )  // dy has changed sign, too
			{	std::cout << "I cannot check the orientation if the curve "
			            << "has too sharp angles - sorry" << std::endl;
				// we can do better here
				exit ( 1 );                                                  }
			if ( dx_pos == dy_pos ) counter++;
			else counter--;                                                    }
		dx_pos = dx_pos_now;  dy_pos = dy_pos_now;                              }
	assert ( ( counter == 2 ) or ( counter == -2 ) );
	return counter == 2;                                                         }

//-------------------------------------------------------------------------------------------------

void switch_orientation ( Cell cll )

// this is always called from switch_orientation ( Mesh )
// we should deal with segments separately, as well with negative cells !
	
{	Mesh msh = cll.boundary();
	std::vector < Cell > vec_of_cells;
	vec_of_cells.reserve ( msh.number_of ( tag::cells_of_dim, msh.dim() ) );
	CellIterator itt = msh.iter_over ( tag::cells_of_dim, msh.dim(), tag::force_positive );
	for ( itt.reset(); itt.in_range(); itt++ )
		vec_of_cells.push_back ( *itt );
  for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from ( msh );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->reverse().add_to ( msh );                                                                 }

//-------------------------------------------------------------------------------------------------

void switch_orientation ( Mesh msh )

// do not use reverse !
// call switch_orientation on each cell

{	std::vector < Cell > vec_of_cells;
	vec_of_cells.reserve ( msh.number_of ( tag::cells_of_dim, msh.dim() ) );
	CellIterator itt = msh.iter_over ( tag::cells_of_dim, msh.dim() );
	for ( itt.reset(); itt.in_range(); itt++ )
		vec_of_cells.push_back ( *itt );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->remove_from ( msh );
	for ( std::vector<Cell>::iterator it = vec_of_cells.begin(); it != vec_of_cells.end(); it++ )
		it->reverse().add_to ( msh );                                                                 }

//-------------------------------------------------------------------------------------------------

std::vector < double > compute_tangent_vec
( Cell start, double desired_length, bool check_orth, std::vector < double > given_vec )

// computes a vector tangent to Manifold::working at point 'start'

// if third argument is true, candidates will be projected onto the space orthogonal to given_vec
// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length
	
{	// Manifold::Implicit * m_impl =  dynamic_cast<Manifold::Implicit*> ( Manifold::working.core );
	// assert ( m_impl );
	// Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*> ( m_impl->surrounding_space.core );
	// assert ( m_euclid );
	std::vector < double > best_tangent;  double longest_projection = 0.;
	const double len_sq = desired_length * desired_length;
	const size_t nc = Manifold::working.coordinates().nb_of_components();
	for ( size_t n = 1; n <= nc; n++ )
	{	const double coef = desired_length / std::sqrt(n);
		// we make sums of n vectors in the canonical basis with both signs
		std::vector < size_t > indices ( n+1 );
		for ( size_t i = 0; i < n; i++ )  indices[i] = i;
		indices[n] = nc;
		while ( true )
		{	std::vector < short int > signs ( n, 1. );
			while ( true )
			{	std::vector < double > tangent ( nc, 0. );
				for ( size_t i = 0; i < n; i++ )  tangent[indices[i]] = signs[i];
				// we normalize 'tangent'
				for ( size_t i = 0; i < nc; i++ ) tangent[i] *= coef;
				// we project
				for ( size_t i = 0; i < nc; i++ )
				{	Function x = Manifold::working.coordinates()[i];
					x ( Mesh::Progressive::temporary_vertex ) = x ( start ) + tangent[i];  }
				Manifold::working.project ( Mesh::Progressive::temporary_vertex );
				for ( size_t i = 0; i < nc; i++ )
				{	Function x = Manifold::working.coordinates()[i];
					tangent[i] = x ( Mesh::Progressive::temporary_vertex ) - x ( start );  }
				if ( check_orth )
				{	double prod = 0.;
					for ( size_t i = 0; i < nc; i++ )
					{	Function x = Manifold::working.coordinates()[i];
						prod += tangent[i] * given_vec[i];                 }
					double lambd = prod / len_sq;
					for ( size_t i = 0; i < nc; i++ )  tangent[i] -= lambd * given_vec[i];  }
				// we choose the longest projection
				double n2 = 0.;
				for ( size_t i = 0; i < nc; i++ )
				{	double tmp = tangent[i];  n2 += tmp*tmp;  }
				if ( n2 > longest_projection )
				{	best_tangent = tangent;  longest_projection = n2;  }
				// now change signs
				bool found = false;
				for ( int i = n-1; i >= 0; i-- )
					if ( signs[i] == 1 )
					{	found = true;  signs[i] = -1;
						for ( size_t j = i+1; j < n; j++ ) signs[j] = 1;
						break;                                            }
				if ( not found ) break;                                                  }
 			// now change indices
			bool found = false;
			for ( int i = n-1; i >= 0; i-- )
				if ( indices[i] < indices[i+1] - 1 )
				{	found = true;
					indices[i]++;
					for ( size_t j = i+1; j<n; j++ )
					{	indices[j] = indices[j-1] + 1;
						assert ( indices[j] < nc );     }
					break;                                                            }
			if ( not found ) break;                                                  }
	}  // end of  for n
	// normalize best_tangent
	double n2 = 0.;
	for ( size_t i = 0; i < nc; i++ )
	{	double tmp = best_tangent[i];  n2 += tmp*tmp;  }
	double norm = approx_sqrt ( n2, tag::around, len_sq, desired_length );
	double coef = desired_length / approx_sqrt ( n2, tag::around, norm*norm, norm );
	for ( size_t i = 0; i < nc; i++ )  best_tangent[i] *= coef;
	return best_tangent;

}  // end of  compute_tangent_vec

//-------------------------------------------------------------------------------------------------

inline std::vector < double > compute_tangent_vec
( const tag::StartAt &, Cell start, const tag::DesiredLength &, double desired_length,
  const tag::OrthogonalTo &, std::vector < double > given_vec                          )

// computes a vector tangent to Manifold::working at point 'start', normal to given_vec

// given_vec must be tangent to Manifold::working at point 'start'
// and must have length approximately equal to desired_length

{	return compute_tangent_vec ( start, desired_length, true, given_vec );  }
	// 'true' as third argument means "do check orthogonality"
	
//-------------------------------------------------------------------------------------------------

inline std::vector < double > compute_tangent_vec
( const tag::StartAt &, Cell start, const tag::DesiredLength &, double desired_length )

// computes a vector tangent to Manifold::working at point 'start'

{	return compute_tangent_vec ( start, desired_length, false, std::vector<double>() );  }
	// 'false' as third argument means "do not check orthogonality"

//-------------------------------------------------------------------------------------------------

void progressive_construct ( Mesh & msh,
	const tag::StartAt &, const Cell & start, const tag::Towards &, std::vector<double> & tangent,
	const tag::StopAt &, const Cell & stop, const tag::DesiredLength &, double desired_length )

// builds a one-dimensional mesh (a curve)
	
// 'start' and 'stop' are vertices (may be one and the same)
	
{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	double len_sq = desired_length * desired_length;
	double augm_length = desired_length * 1.618034,  // golden number
	       augm_len_sq = augm_length * augm_length;
	size_t counter = 1;
	size_t max_counter = 0;  // std::cin >> max_counter;
	Cell A = start;
	while ( true )
	{	double d = Manifold::working.dist_sq ( A, stop );
		if ( d < augm_len_sq )
		{	std::vector < double > e ( nc );
			double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				e[i] = x(stop) - x(A);
				prod += tangent[i] * e[i];                        }
			if ( prod > 0. )
			{	Cell last ( tag::segment, A.reverse(), stop );
				last.add_to ( msh );
				// redistribute vertices
				double n2 = 0.;
				for ( size_t i = 0; i < nc; i++ )
				{	double temp = e[i];  temp *= temp;  n2 += temp;  }
				double norm = approx_sqrt ( n2, tag::around, len_sq, desired_length );
				norm = approx_sqrt ( n2, tag::around, norm*norm, norm );
				redistribute_vertices ( msh, start, stop, desired_length, norm, 15 );
				break;                                                                 }  }
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x(B) = x(A) + tangent[i];                         }
		if ( counter == max_counter ) return;
		Manifold::working.project ( B );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tangent[i] = x(B) - x(A);                         }
		double n2 = 0.;
		for ( size_t i = 0; i < nc; i++ )
		{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
		n2 = approx_sqrt ( n2, tag::around, len_sq, desired_length );
		n2 = desired_length / n2;
		for ( size_t i = 0; i < nc; i++ )  tangent[i] *= n2;
		Cell AB ( tag::segment, A.reverse(), B );
		AB.add_to ( msh );  counter++;  A = B;                                 }

} // end of  progressive_construct

//-------------------------------------------------------------------------------------------------

void progressive_construct ( Mesh & msh, const tag::StartAt &, const Cell & start,
	const tag::StopAt &, const Cell & stop, const tag::DesiredLength &, double desired_length )
	
// builds a one-dimensional mesh (a curve)
	
// no tangent vector provided, thus no way of knowing which way to go
// the only solution is to start walking in both directions simultaneously
// the first one to reach 'stop' wins

// 'start' and 'stop' are vertices (may be one and the same)

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	std::vector < double > best_tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, desired_length );

	// start walking along the manifold from 'start' in the direction of best_tangent
	// and, simultaneously, in the opposite direction, given by -best_tangent
	double augm_length = desired_length * 1.5,
	       augm_len_sq = augm_length * augm_length;
	std::vector < double > tan1 = best_tangent, tan2 = best_tangent;
	for ( size_t i = 0; i < nc; i++ ) tan2[i] *= -1.;
	Cell ver1 ( tag::vertex ), ver2 ( tag::vertex );
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x ( ver1 ) = x ( start );  x ( ver2 ) = x ( start );  }
	int winner;  //  will be 1 or -1
	while ( true )
	{	for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( Mesh::Progressive::temporary_vertex ) = x ( ver1 ) + tan1[i];  }
		Manifold::working.project ( Mesh::Progressive::temporary_vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tan1[i] = x ( Mesh::Progressive::temporary_vertex ) - x ( ver1 );
			x ( ver1 ) = x ( Mesh::Progressive::temporary_vertex );            }
		double d = Manifold::working.dist_sq ( ver1, stop );
		if ( d < augm_len_sq )
		{	double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				prod += tan1[i] * ( x(stop) - x(ver1) );          }
			if ( prod > 0. )  { winner = 1;  break;  }           }
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( Mesh::Progressive::temporary_vertex ) = x ( ver2 ) + tan2[i];  }
		Manifold::working.project ( Mesh::Progressive::temporary_vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tan2[i] = x ( Mesh::Progressive::temporary_vertex ) - x ( ver2 );
			x ( ver2 ) = x ( Mesh::Progressive::temporary_vertex );            }
		d = Manifold::working.dist_sq ( ver2, stop );
		if ( d < augm_len_sq )
		{	double prod = 0.;
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				prod += tan2[i] * ( x(stop) - x(ver2) );          }
			if ( prod > 0. )  { winner = -1;  break;  }           }
	}  // end of  while true

	ver1.dispose();  ver2.dispose();

	for ( size_t i = 0; i < nc; i++ ) best_tangent[i] *= winner;
	progressive_construct ( msh, tag::start_at, start, tag::towards, best_tangent,
	                        tag::stop_at, stop, tag::desired_length, desired_length   );
}

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start, const tag::StopAt &, const Cell & stop,
  const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & )

// 'start' and 'stop' are vertices (may be one and the same)

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	assert ( start.dim() == 0 );
	assert ( stop.dim() == 0 );

	std::vector < double > best_tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, desired_length );

	Mesh msh1 ( tag::of_dimension_one );
	progressive_construct ( msh1, tag::start_at, start, tag::towards, best_tangent,
                          tag::stop_at, stop, tag::desired_length, desired_length );
	for ( size_t i = 0; i < nc; i++ )  best_tangent[i] *= -1.;
	Mesh msh2 ( tag::of_dimension_one );
	progressive_construct ( msh2, tag::start_at, start, tag::towards, best_tangent,
                          tag::stop_at, stop, tag::desired_length, desired_length );
	switch_orientation ( msh2 );
	Mesh whole ( tag::join, msh1, msh2 );

	if ( correctly_oriented ( whole ) )
	{	msh2.dispose();  msh = msh1;  }
	else
	{	msh1.dispose();  switch_orientation ( msh2 );  msh = msh2;  }

	whole.dispose();
}
	
//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::StartAt &, Cell start,
  const tag::Towards &, std::vector<double> & normal,
  const tag::Boundary &, Mesh bdry,
  const tag::DesiredLength, double desired_length        )

// for two-dimensional meshes (arbitrary geometric dimension)
	
// 'start' is a vertex or segment belonging to 'bdry'
// 'normal' is a vector tangent to the working manifold, orthogonal to 'start'

{	Mesh::Progressive::desired_length = desired_length;
	Mesh::Progressive::desired_length_sq = desired_length * desired_length;
	Mesh::Progressive::short_dist =
		Mesh::Progressive::desired_length / Mesh::Progressive::one_plus_tolerance;
	Mesh::Progressive::short_dist_sq =
		Mesh::Progressive::short_dist * Mesh::Progressive::short_dist;
	Mesh::Progressive::long_dist =
		Mesh::Progressive::desired_length * Mesh::Progressive::one_plus_tolerance;
	Mesh::Progressive::long_dist_sq =
		Mesh::Progressive::long_dist * Mesh::Progressive::long_dist;

	// we don't want to change 'bdry' so we make a copy of it
	// Mesh interface ( tag::deep_copy_of, bdry );  // wait, we want to switch to a list of interfaces ...
	Mesh interface = bdry;
	Mesh::Progressive::interface = interface;
	Mesh::Progressive::under_constr = msh;
	const size_t nc = Manifold::working.coordinates().nb_of_components();
	Mesh::Progressive::nb_of_coords = nc;
	Cell vertex_recently_built ( tag::non_existent );
	std::set < Cell::Core* > set_of_nearby_vertices;
	// vertices close to vertex_recently_built

	Mesh::Cloud cl ( tag::min_dist, desired_length, tag::ratio, 6. );

	{ // just a block of code for hiding variables
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[0];
	std::list<Cell::Core*>::iterator it = l.begin();
	for ( ; it != l.end(); it++ )
	{	Cell::Core * ver_p = *it;
		Cell ver ( tag::whose_core_is, ver_p );
		Mesh::Cloud::add_point ( ver );          }
	} // just a block of code for hiding variables

	if ( start.dim() != 1 )
	{	assert ( start.dim() == 0 );
		for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			std::cout << x(start) << " ";                    }
		std::cout << std::endl;
		start = interface.cell_in_front_of ( start );                    }
	assert ( start.dim() == 1 );
	assert ( bdry.dim() == 1 );
	assert ( start.belongs_to ( interface, tag::oriented ) );
	
	{ // just a block of code for hiding variables
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = normal[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  normal[i] *= coef;
	Cell start_base = start.base().reverse();
	Cell start_tip = start.tip();
	std::vector < double > vec ( nc );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		vec[i] = x ( start_tip ) - x ( start_base );      }
	Mesh::Progressive::improve_normal ( start_base, start_tip, vec, normal );
	// ensures again the norm is right, projects on the tangent space,
	// ensures orthogonality with 'start' segment
	} // just a block of code for hiding variables

	Mesh::Progressive::normals[start.core] = normal;
	Cell ret = Mesh::Progressive::build_normals ( start );
	assert ( ret == start );

	Cell point_60 = start.tip();
	
	int stopping_criterion = 0;  // std::cin >> stopping_criterion;
	int current_name = 1;

restart:
	
	Cell stop_point_60 = point_60, point_120 = point_60, stop_point_120 = point_60;

angles_60 :

	// starting at 'point_60', we go along this connected component of 'interface'
	{ // just a block of code for hiding variables
	Cell prev_seg = interface.cell_behind ( point_60, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	Cell next_seg = interface.cell_in_front_of ( point_60, tag::surely_exists );
		Cell B = next_seg.tip();  assert ( A != B );
		double d = Manifold::working.dist_sq ( A, B );
		// assert ( d > Mesh::Progressive::short_dist_sq );
		if ( d < Mesh::Progressive::long_dist_sq )  // we may want to form a triangle
		// but first we must make sure it is on the correct side
		if ( Mesh::Progressive::positive_orientation ( A, point_60, prev_seg, next_seg ) )
		if ( Mesh::Progressive::cos_sq_60 ( A, point_60, B, prev_seg, next_seg) > 0.03 )
		// triangle waiting to be filled; see paragraph 9.7 in the manual
		{	Cell seg_next_to_B = interface.cell_in_front_of(B);
			Cell ver_next_to_B = seg_next_to_B.tip();
			set_of_nearby_vertices.erase ( point_60.core );
			if ( ver_next_to_B == A )  // this is the last triangle in this piece of interface
			{	Mesh::Progressive::fill_last_triangle
					( A, point_60, B, prev_seg, next_seg, seg_next_to_B, interface );
				std::cout << "shrinking triangle " << ++current_name << std::endl << std::flush;
				if ( current_name == stopping_criterion ) return;			
				goto search_for_start;  	                                                        }
			Cell AB ( tag::segment, A.reverse(), B );
			Mesh::Progressive::fill_60 ( prev_seg, next_seg, AB.reverse(), point_60, interface );
			AB.add_to ( interface );
			Mesh::Progressive::build_one_normal ( A, B, AB );  // based on previous segment
			std::cout << "found angle around 60 deg " << ++current_name << std::endl << std::flush;
			if ( current_name == stopping_criterion ) return;
			if ( stop_point_120 == B )
			{	if ( stop_point_120 == point_120 )
				{	point_120 = A;  stop_point_120 = A;  }
				else  stop_point_120 = ver_next_to_B;    }
			if ( stop_point_120 == point_60 )  // we have all the loop to cover for 120 deg
			{	assert ( point_120 == stop_point_120 );
				stop_point_120 = A;                     }
			if ( point_120 == point_60 ) point_120 = A;
			if ( stop_point_60 == B )  stop_point_60 = ver_next_to_B;
			if ( stop_point_60 == point_60 )  // we have all the loop to cover for 60 deg
				stop_point_60 = A;
			point_60 = A;  goto angles_60;                                                    }
		if ( B == stop_point_60 ) break;
		A = point_60;  prev_seg = next_seg;  point_60 = B;                                            }
	} // just a block of code for hiding variables

check_touching :

	{ // just a block of code for hiding 'touch'
	bool touch = Mesh::Progressive::check_touching
		( vertex_recently_built, set_of_nearby_vertices );
	vertex_recently_built = Cell ( tag::non_existent );
	if ( touch )
	{	assert ( point_120.belongs_to ( interface, tag::not_oriented ) );
		std::cout << "touch " << ++current_name << std::endl << std::flush;
		if ( current_name == stopping_criterion ) return;			
		point_60 = point_120;  stop_point_60 = point_60;  stop_point_120 = point_60;
		goto angles_60;                                                               }
	} // just a block of code for hiding 'touch'

// look for angles around 120 deg :

	// starting at 'point_120', we go along this connected component of 'interface'
	{  // just a block of code for hiding prev_seg and A
	Cell prev_seg = interface.cell_behind ( point_120, tag::surely_exists );
	Cell A = prev_seg.base().reverse();
	while ( true )
	{	Cell next_seg = interface.cell_in_front_of ( point_120, tag::surely_exists );
		Cell  B = next_seg.tip();
		if ( Mesh::Progressive::cos_sq_120 ( A, point_120, B, prev_seg, next_seg) < 0.55 )  // 0.67
		// angle around 120 deg, we want to form two triangles; see paragraph 9.7 in the manual
		{	// we don't build a new vertex yet, we want to check for a quadrangle first
			Cell seg_prev_to_A = interface.cell_behind ( A );
			Cell seg_next_to_B = interface.cell_in_front_of ( B );
			Cell ver_prev_to_A = seg_prev_to_A.base().reverse();
			Cell ver_next_to_B = seg_next_to_B.tip();
			if ( ver_prev_to_A == ver_next_to_B )  // quadrangle
			{	// check orientations : correct side ?
				std::cout << "shrinking quadrangle" << std::endl << std::flush;
				// choose the shortest diagonal, add two triangles to the mesh under construction
				if ( Manifold::working.dist_sq ( A, B ) <
						 Manifold::working.dist_sq ( point_120, ver_next_to_B ) )
				{	Cell AB ( tag::segment, A.reverse(), B );
					Mesh::Progressive::fill_60 ( prev_seg, next_seg, AB.reverse(), point_120, interface );
					Mesh::Progressive::fill_60 ( seg_next_to_B, seg_prev_to_A, AB, ver_prev_to_A, interface ); }
				else
				{	Cell seg ( tag::segment, point_120.reverse(), ver_next_to_B );
					Mesh::Progressive::fill_60 ( next_seg, seg_next_to_B, seg.reverse(), B, interface );
					Mesh::Progressive::fill_60 ( seg_prev_to_A, prev_seg, seg, A, interface );            }
				goto search_for_start;                                                                          }
			Cell P ( tag::vertex );  vertex_recently_built = P;
			// now we want to place this new vertex accordingly
			assert ( Mesh::Progressive::normals.find(prev_seg.core) != Mesh::Progressive::normals.end() );
			assert ( Mesh::Progressive::normals.find(next_seg.core) != Mesh::Progressive::normals.end() );
			std::vector < double > & nor_a = Mesh::Progressive::normals [ prev_seg.core ];
			std::vector < double > & nor_b = Mesh::Progressive::normals [ next_seg.core ];
			std::vector < double > sum_of_nor ( nc );
			for ( size_t i = 0; i < nc; i++ )  sum_of_nor[i] = nor_a[i] + nor_b[i];
			for ( size_t i = 0; i < nc; i++ )
			{	Function x = Manifold::working.coordinates()[i];
				x(P) = x(point_120) / 2. + ( x(A) + x(B) ) / 4.
					+ sum_of_nor[i] * Mesh::Progressive::half_of_sqrt_of_075;  }
			Manifold::working.project(P);
			Cell AP ( tag::segment, A.reverse(), P );
			Cell BP ( tag::segment, B.reverse(), P );
			Cell PB = BP.reverse();
			Cell sP ( tag::segment, point_120.reverse(), P );
			Cell tri1 ( tag::triangle, prev_seg, sP, AP.reverse() );
			Cell tri2 ( tag::triangle, next_seg, BP, sP.reverse() );
			tri1.add_to ( msh );
			tri2.add_to ( msh );
			prev_seg.remove_from ( interface );
			next_seg.remove_from ( interface );
			Mesh::Progressive::normals.erase ( prev_seg.core );
			Mesh::Progressive::normals.erase ( next_seg.core );
			AP.add_to ( interface );
			PB.add_to ( interface );
			Mesh::Cloud::remove_point ( point_120 );
			Mesh::Progressive::build_one_normal ( A, P, AP );  // based on previous segment
			Mesh::Progressive::build_one_normal ( P, B, PB );  // based on previous segment
			Mesh::Progressive::relocate ( P, 2, sum_of_nor, nc, set_of_nearby_vertices );
			// find more vertices close to P and take them all into account; modifies sum_of_nor
			assert ( prev_seg.tip() == point_120 );
		  msh.baricenter ( point_120, prev_seg );
			Mesh::Cloud::add_point ( P );
			std::cout << "found angle around 120 deg " << ++current_name << std::endl << std::flush;
			if ( current_name == stopping_criterion ) return;			
			if ( stop_point_120 == point_120 )  // we have all the loop to cover
				stop_point_120 = A;
			if ( stop_point_120 == B )  stop_point_120 = ver_next_to_B;
			point_120 = A;
			point_60 = A;
			stop_point_60 = ver_next_to_B;
			goto angles_60;                                                        }
		if ( B == stop_point_120 ) break;
		A = point_120;  prev_seg = next_seg;  point_120 = B;
	}  // end of while
	}  // just a block of code for hiding prev_seg and A

// build a brand new triangle :
	
	{ // just a block of code for hiding variables
	// we use point_120, any other point would do
	Cell next_seg = interface.cell_in_front_of ( point_120, tag::surely_exists );
	Cell B = next_seg.tip();
	Cell P ( tag::vertex );  vertex_recently_built = P;
	assert ( Mesh::Progressive::normals.find(next_seg.core) != Mesh::Progressive::normals.end() );
	std::vector < double > & f = Mesh::Progressive::normals [ next_seg.core ];
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = ( x(point_120) + x(B) ) / 2. + f[i] * Mesh::Progressive::sqrt_of_075;  }
	Manifold::working.project(P);
	Cell AP ( tag::segment, point_120.reverse(), P );
	Cell BP ( tag::segment, B.reverse(), P );
	Cell PB = BP.reverse();
	Cell tri ( tag::triangle, next_seg, BP, AP.reverse() );
	tri.add_to ( msh );
	next_seg.remove_from ( interface );
	Mesh::Progressive::normals.erase ( next_seg.core );
	AP.add_to ( interface );
	PB.add_to ( interface );
	// if ( current_name+1 == stopping_criterion ) return;			
	std::cout << "building brand new triangle " << ++current_name << std::endl << std::flush;
	if ( current_name == stopping_criterion ) return;			
	Mesh::Progressive::build_one_normal ( point_120, P, AP );  // based on previous segment
	Mesh::Progressive::build_one_normal ( P, B, PB );  // based on previous segment
	Mesh::Progressive::relocate ( P, 1, f, nc, set_of_nearby_vertices );
	// find more vertices close to P and take them all into account; modifies f
	Mesh::Cloud::add_point ( P );
	stop_point_120 = interface.cell_in_front_of(B).tip();
	goto check_touching;
	} // just a block of code for hiding variables

search_for_start :
// execution only reaches this point through 'goto'

	{ // just a block of code for hiding variables
	// we look for a segment in 'interface' which has a normal
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[1];
	if ( l.size() == 0 ) return;  // empty interface, meshing process ended
	std::cout << "search for start " << l.size() << std::endl << std::flush;
	Cell::Core * start_seg_p = nullptr;
	assert ( l.size() >= 2 );
	std::list<Cell::Core*>::iterator it = l.begin();
	for ( ; it != l.end(); it++ )
		if ( Mesh::Progressive::normals.find ( *it ) != Mesh::Progressive::normals.end() )
		{	start_seg_p = *it;  break;  }
	assert ( start_seg_p );
	point_60.core = start_seg_p->tip();
	// any point on this connected component would do
	goto restart;
	} // just a block of code for hiding variables

}  // end of progressive_construct

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::DesiredLength &, double desired_length, bool check_and_switch )
	
// last argument tells whether to check the orientation of the resulting mesh
// and switch it if necessary

{	const size_t nc = Manifold::working.coordinates().nb_of_components();
	
	// call to 'search_start_ver' does not depend on the dimension of the mesh
	Cell start = search_start_ver ( desired_length );

	// call to 'compute_tangent_vec' does not depend on the dimension of the mesh
	std::vector < double > tangent = compute_tangent_vec
		( tag::start_at, start, tag::desired_length, desired_length );

	// now we branch, depending on the dimension
	if ( msh.dim() == 1 )
		progressive_construct ( msh, tag::start_at, start, tag::towards, tangent,
		                        tag::stop_at, start, tag::desired_length, desired_length );
	else
	{	assert ( msh.dim() == 2 );  // no 3D for now
		assert ( nc == 3 );
		Cell B ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( B ) = x ( start ) + tangent[i];              }
		Cell AB ( tag::segment, start.reverse(), B );
		std::vector < double > normal = compute_tangent_vec ( tag::start_at, start,
		  tag::desired_length, desired_length, tag::orthogonal_to, tangent );
		Cell C ( tag::vertex );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			x ( C ) = x ( start ) + 0.5 * tangent[i] - 0.866 * normal[i];  }
		Cell BC ( tag::segment, B.reverse(), C );
		Cell CA ( tag::segment, C.reverse(), start );
		Cell tri ( tag::triangle, AB, BC, CA );
		tri.add_to ( msh );
		Mesh interf ( tag::of_dimension_one );
		AB.reverse().add_to(interf);
		BC.reverse().add_to(interf);
		CA.reverse().add_to(interf);
		progressive_construct ( msh, tag::start_at, AB.reverse(), tag::towards, normal,
		                        tag::boundary, interf, tag::desired_length, desired_length );
	}  // end of  else  with  msh.dim() == 2

	if ( not check_and_switch ) return;
	if ( not correctly_oriented ( msh ) )  switch_orientation ( msh );
	assert ( correctly_oriented ( msh ) );                                                }

//-------------------------------------------------------------------------------------------------
	
void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double desired_length,
	const tag::InherentOrientation, bool check_and_switch )

{	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	assert ( msh.dim() == 2 );
	assert ( start.core->belongs_to ( interface.core, tag::oriented ) );
	Mesh interf_rev ( tag::of_dimension_one );
	CellIterator it = interface.iter_over ( tag::cells_of_dim, 1 );
	for ( it.reset(); it.in_range(); it++ )
		(*it).reverse().add_to ( interf_rev );
	assert ( start.reverse().core->belongs_to ( interf_rev.core, tag::oriented ) );

	// compute a normal vector, on an arbitrary side of 'start'
	Cell A = start.base().reverse();
	Cell B = start.tip();
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	Function z = Manifold::working.coordinates()[2];
	std::vector < double > tan { x(B) - x(A), y(B) - y(A), z(B) - z(A) };
	std::vector < double > nor = compute_tangent_vec
		( tag::start_at, A, tag::desired_length, desired_length, tag::orthogonal_to, tan );

	progressive_construct ( msh, tag::start_at, start, tag::towards, nor,
	                        tag::boundary, interface, tag::desired_length, desired_length );	

	if ( not check_and_switch ) return;

	// build the mesh on the other side of 'interface'
	for ( size_t i = 0; i < 3; i++ )  nor[i] *= -1.;
	Mesh msh2 ( tag::of_dimension, 2, tag::greater_than_one );
	progressive_construct ( msh2, tag::start_at, start.reverse(), tag::towards, nor,
	                        tag::boundary, interf_rev, tag::desired_length, desired_length );	
	// join everything to get a mesh on the entire manifold
	std::cout << "wait a minute ..." << std::endl;
	// we should switch_orientation ( msh2 )
	Mesh glob ( tag::join, msh, msh2 );

	if ( not correctly_oriented ( glob ) )
	{	switch_orientation ( msh2 );  msh = msh2;  msh.dispose();  }
	else  msh2.dispose();
	glob.dispose();
}

//-------------------------------------------------------------------------------------------------
	
void progressive_construct
( Mesh & msh, const tag::StartAt &, const Cell & start,
  const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double desired_length        )

// for two-dimensional meshes in RR^2 (intrinsic orientation)
//   or in a 2D submanifold of RR^3 (random orientation)
	
// 'start' is a segment belonging to 'interface'
// no normal vector provided, we need to build our own

{	if ( Manifold::working.coordinates().nb_of_components() == 3 )
	// surface in RR^3
	{	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
		        tag::desired_length, desired_length, tag::inherent_orientation, false );
		// last argument 'false' means do not check orientation, leave it random
		return;                                                                            }
	
	// domain in the plane RR^2
	Cell A = start.base().reverse();
	Cell B = start.tip();
	Function x = Manifold::working.coordinates()[0];
	Function y = Manifold::working.coordinates()[1];
	std::vector < double > nor { y(A) - y(B), x(B) - x(A) };   // rotate with 90 deg

	progressive_construct ( msh, tag::start_at, start, tag::towards, nor,
	                        tag::boundary, interface, tag::desired_length, desired_length );	
}
	
//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double desired_length    )

// for two-dimensional meshes in RR^2 (intrinsic orientation)
//   or in a 2D submanifold of RR^3 (random orientation)
	
{	// we search for a starting point
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[1];
	std::list<Cell::Core*>::iterator it = l.begin();
	Cell start ( tag::whose_core_is, *it );

	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
	                        tag::desired_length, desired_length                  );
}

//-------------------------------------------------------------------------------------------------

void progressive_construct
( Mesh & msh, const tag::Boundary &, Mesh interface,
  const tag::DesiredLength, double desired_length,
	const tag::InherentOrientation, bool check_and_switch )

// if last argument is true, compute inherent orientation
// otherwise, random orientation

{	// we search for a starting point
	// 'interface' is a one-dimensional mesh, not necessarily connected
	// so we cannot use a CellIterator - perhaps an unstructured one ?
	std::list < Cell::Core* > & l = interface.core->cells[1];
	std::list<Cell::Core*>::iterator it = l.begin();
	Cell start ( tag::whose_core_is, *it );

	progressive_construct ( msh, tag::start_at, start, tag::boundary, interface,
	                        tag::desired_length, desired_length,
                          tag::inherent_orientation, check_and_switch );
}

//-------------------------------------------------------------------------------------------------

inline size_t get_topological_dim ( )
{	Manifold::Implicit::OneEquation * m_impl_1 =
		dynamic_cast<Manifold::Implicit::OneEquation*> ( Manifold::working.core );
	if ( m_impl_1 )
	{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*>
			( m_impl_1->surrounding_space.core );
		assert ( m_euclid );
		return m_euclid->coord_func.nb_of_components() - 1;              }
	Manifold::Implicit::TwoEquations * m_impl_2 =
		dynamic_cast<Manifold::Implicit::TwoEquations*> ( Manifold::working.core );
	if ( m_impl_2 )
	{	Manifold::Euclid * m_euclid = dynamic_cast<Manifold::Euclid*>
			( m_impl_2->surrounding_space.core );
		assert ( m_euclid );
		return m_euclid->coord_func.nb_of_components() - 2;              }
	assert ( false );                                                               } 

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, true );
	// last argument true means : check orientation, switch it if necessary

	Mesh::Progressive::temporary_vertex.dispose();                                }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double desired_length                     )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, true );
	// last argument true means : check orientation, switch it if necessary

	Mesh::Progressive::temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                                 }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length,
             const tag::RandomOrientation &                                               )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, false );
	// last argument false means : do not check orientation, leave it random

	Mesh::Progressive::temporary_vertex.dispose();                                 }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double desired_length, const tag::RandomOrientation & )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	
	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, false );
	// last argument false means : do not check orientation, leave it random

	Mesh::Progressive::temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                                  }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::DesiredLength &, double desired_length,
             const tag::InherentOrientation &                                             )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	assert ( Manifold::working.coordinates().nb_of_components() == this->dim() + 1 );
	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, true );
	// last argument true means : check orientation, switch it if necessary

	Mesh::Progressive::temporary_vertex.dispose();                                      }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
             const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & )

:	Mesh ( tag::of_dimension, get_topological_dim(), tag::might_be_one )

{	assert ( Manifold::working.coordinates().nb_of_components() == this->dim() + 1 );
	Manifold tmp_manif = Manifold::working;
	Manifold::working = manif;
	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this,  tag::desired_length, desired_length, true );
	// last argument true means : check orientation, switch it if necessary

	Mesh::Progressive::temporary_vertex.dispose();
	Manifold::working = tmp_manif;                                                        }

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & tangent,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double desired_length              )

// 'start' and 'stop' may be the same cell

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	const size_t nc = Manifold::working.coordinates().nb_of_components();	
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  tangent[i] *= coef;

	improve_normal ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
                          tag::stop_at, stop, tag::desired_length, desired_length );

	Mesh::Progressive::temporary_vertex.dispose();
}	

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & tangent,
             const tag::DesiredLength &, double desired_length              )

// 'start' and 'stop' may be the same cell

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	const size_t nc = Manifold::working.coordinates().nb_of_components();	
	double n2 = 0.;  // ensure that normal has the right norm
	for ( size_t i = 0; i < nc; i++ )
	{	double temp = tangent[i];  temp *= temp;  n2 += temp;  }
	double coef = desired_length / std::sqrt(n2);
	for ( size_t i = 0; i < nc; i++ )  tangent[i] *= coef;
	
	improve_normal ( start, tangent );
	// ensures again the norm is right, then projects on the tangent space

	progressive_construct ( *this, tag::start_at, start, tag::towards, tangent,
                          tag::stop_at, start, tag::desired_length, desired_length );

	Mesh::Progressive::temporary_vertex.dispose();
}	

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double desired_length )

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,
	                        tag::desired_length, desired_length              );

	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::DesiredLength &, double desired_length                   )

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, start,
	                        tag::desired_length, desired_length               );

	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
             const tag::StopAt &, const Cell & stop,
             const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & )

:	Mesh ( tag::of_dimension_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );
	
	progressive_construct ( *this, tag::start_at, start, tag::stop_at, stop,
	                        tag::desired_length, desired_length, tag::inherent_orientation );

	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double desired_length                 )

// for now, only works for two-dimensional meshes (either in RR2 or in RR3)
// should be adapted for three-dimensional meshes

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this, tag::boundary, interface, tag::desired_length, desired_length );
	
	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double desired_length, const tag::IntrinsicOrientation & )

// for now, only works for two-dimensional meshes in RR2
// should be adapted for three-dimensional meshes

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	assert ( Manifold::working.coordinates().nb_of_components() == 2 );
	progressive_construct ( *this, tag::boundary, interface, tag::desired_length, desired_length );
	
	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::DesiredLength &, double desired_length, const tag::InherentOrientation & )

// for two-dimensional meshes in RR^3

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	assert ( Manifold::working.coordinates().nb_of_components() == 3 );
	progressive_construct ( *this, tag::boundary, interface,
                          tag::desired_length, desired_length, tag::inherent_orientation, true );
	
	Mesh::Progressive::temporary_vertex.dispose();
}

//-------------------------------------------------------------------------------------------------

Mesh::Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
             const tag::StartAt &, const Cell & start,
             const tag::Towards &, std::vector<double> & normal,             
             const tag::DesiredLength &, double desired_length                )
	
// 'start' is a vertex or segment belonging to 'interface'

:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )  // positive, by default

{	Mesh::Progressive::temporary_vertex = Cell ( tag::vertex );

	progressive_construct ( *this, tag::start_at, start, tag::towards, normal,
													tag::boundary, interface, tag::desired_length, desired_length );
	
	Mesh::Progressive::temporary_vertex.dispose();
}
			  
//-------------------------------------------------------------------------------------------------

inline bool Mesh::Progressive::positive_orientation  // static
( const Cell & A, const Cell & B, const Cell & AB, const Cell & BC )
	
{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	// code below is identical to part of Mesh::Progressive::cos_sq_120
	// if you change anything, please change both; keep them identical
	assert ( Mesh::Progressive::normals.find(BC.core) != Mesh::Progressive::normals.end() );
	std::vector < double > e ( Mesh::Progressive::nb_of_coords ),
		& f = Mesh::Progressive::normals[BC.core];
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );              }
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of Mesh::Progressive::cos_sq_120
	return prod < 0.;                                                                      }

//-------------------------------------------------------------------------------------------------

inline double Mesh::Progressive::cos_sq_60  // static
( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )

// return cosine square of 180 - ABC (or 0. if wrong orientation)
// check that  cos_sq > 0.03  to ensure angle ABC is below 80 deg

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	std::vector < double > e1 ( Mesh::Progressive::nb_of_coords ),
	                       e2 ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e1[i] = x ( B ) - x ( A );
		e2[i] = x ( C ) - x ( B );                        }
	double prod = 0.;  // scalar product e1.e2
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  prod += e1[i] * e2[i];
	// use Riemannian product !
	if ( prod > 0. ) return 0.;
	// use Riemannian metric below
	double norm1 = 0., norm2 = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	norm1 += e1[i] * e1[i];
		norm2 += e2[i] * e2[i];  }
	return prod * prod / norm1 / norm2;  // cosine square
}

//-------------------------------------------------------------------------------------------------
			
inline double Mesh::Progressive::cos_sq_120  // static
( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & BC )

// return cosine square of 180 - ABC (or 2. if wrong orientation)
// for instance, calling function should check cos_sq < 0.671 to ensure angle ABC is below 145 deg
// or check cos_sq < -0.17365 to ensure angle ABC is below 80 deg

{	assert ( A == AB.base().reverse() );
	assert ( B == AB.tip() );
	assert ( B == BC.base().reverse() );
	assert ( C == BC.tip() );
	// first, check orientation
	// code below is identical to part of Mesh::Progressive::positive_orientation
	// if you change anything, please change both; keep them identical
	assert ( Mesh::Progressive::normals.find(BC.core) != Mesh::Progressive::normals.end() );
	std::vector < double > e ( Mesh::Progressive::nb_of_coords ),
		& f = Mesh::Progressive::normals[BC.core];
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e[i] = x ( B ) - x ( A );                         }
	double prod = 0.;  // scalar product e.f
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  prod += e[i] * f[i];
	// we could have used the Riemannian product, but the sign should be the same
	// code above is identical to part of Mesh::Progressive::positive_orientation
	if ( prod > 0. ) return 2.;
	std::vector < double > e2 ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		e2[i] = x ( C ) - x ( B );                       }
	// use Riemannian metric below
	double norm1 = 0., norm2 = 0.;  prod = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
	{	norm1 += e [i] * e [i];
		norm2 += e2[i] * e2[i];
		prod  += e [i] * e2[i];  }
	return prod * prod / norm1 / norm2;  // cosine square
}

//-------------------------------------------------------------------------------------------------
			
void Mesh::Progressive::relocate  // static
( const Cell & P, size_t n, std::vector<double> & normal_dir,
  const size_t & nc, std::set<Cell::Core*> & set_of_ver       )

// re-compute the placement of a newly created vertex

// vertex has been located according to two segments, from angles_120 :   n == 2
// or according to only one segment, if built from a brand new triangle : n == 1
	
// we compute here 'set_of_ver' (which is the set of all vertices in the cloud
// close enough to 'ver') and keep it for future use in 'check_touching'

{	Mesh interface = Mesh::Progressive::interface;

	// make a list (using a set) of nearby points
	// build a vector of segments from it
	// relocate point P by averaging all normals

	set_of_ver = Mesh::Cloud::find_close_vertices ( P, Mesh::Progressive::desired_length );
	// P has not been added to the cloud yet, so it will not show up in 'set_of_ver'
	std::vector < Cell::Core* > vector_of_seg;
	// vector_of_seg will contain all segments whose both extremities belong to 'set_of-ver'
	// but not the segments adjacent to P (since P does not belong to 'set_of_ver')
	for ( std::set<Cell::Core*>::iterator it = set_of_ver.begin(); it != set_of_ver.end(); it++ )
	{	Cell A ( tag::whose_core_is, *it );
		Cell AB = interface.cell_in_front_of ( A );
		if ( set_of_ver.find ( AB.tip().core ) != set_of_ver.end() )
		{	vector_of_seg.push_back ( AB.core );                               }  }
	
	if ( vector_of_seg.empty() ) return;

	// there are two cases : 
	// the piece of the interface that we just encountered may have normals or not
	size_t counter = 0;
	Cell kept_seg ( tag::non_existent );
	for ( size_t i = 0; i < vector_of_seg.size(); i++ )
	{	Cell::Core * seg_p = vector_of_seg[i];
		if ( Mesh::Progressive::normals.find ( seg_p ) == Mesh::Progressive::normals.end() )
		{	counter++; kept_seg.core = seg_p;  }                                             }
	if ( counter > 0 )  // there are 'counter' segments with no normal
	{	assert ( counter == 1 );
		// build normal of 'kept_seg' from 'normal_dir'
		Cell A = kept_seg.base().reverse();
		Cell B = kept_seg.tip();
		std::vector < double > tangent_dir ( nc );
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			tangent_dir[i] = x(B) - x(A);
			normal_dir[i] *= -1.;                             }
		Mesh::Progressive::improve_normal ( A, B, tangent_dir, normal_dir );  // modifies normal_dir
		Mesh::Progressive::normals[kept_seg.core] = normal_dir;
		Cell ret = Mesh::Progressive::build_normals ( kept_seg );
		assert ( ret == kept_seg );                                           }

	std::vector < double > pos = Manifold::working.coordinates() ( P );
	if ( n != 1 )
	{	assert ( n == 2 );
		for ( size_t i = 0; i < nc; i++ )  pos[i] *= n;  }
	for ( size_t j = 0; j < vector_of_seg.size(); j++ )
	{	Cell AB ( tag::whose_core_is, vector_of_seg[j] );
	  std::map < Cell::Core*, std::vector<double> > :: iterator it =
			Mesh::Progressive::normals.find ( AB.core );
		assert ( it != Mesh::Progressive::normals.end() );
		std::vector < double > & nor = it->second;
		// std::vector < double > & nor = Mesh::Progressive::normals [ AB.core ];
		Cell A = AB.base().reverse();
		Cell B = AB.tip();
		for ( size_t i = 0; i < nc; i++ )
		{	Function x = Manifold::working.coordinates()[i];
			pos[i] += ( x(A) + x(B) ) / 2. + nor[i] * Mesh::Progressive::sqrt_of_075;  }  }
	n += vector_of_seg.size();
	// std::cout << "relocate " << n << std::endl << std::flush;
	for ( size_t i = 0; i < nc; i++ )
	{	Function x = Manifold::working.coordinates()[i];
		x(P) = pos[i] / n;                                 }
	
	Manifold::working.project(P);

	// check P is on the right side of kept_seg
}

//-------------------------------------------------------------------------------------------------

inline bool Mesh::Progressive::check_touching  // static
( Cell & ver, std::set<Cell::Core*> & set_of_ver )

// analyse position of recently built vertex 'ver' relatively to other vertices on the interface
// occasionally, different connected components of the interface touch and merge
// or, the current connected component of the interface may touch itself and split in two

// return true if a touch was detected and the corresponding pieces have been merged

// we take advantage of 'set_of_ver' which is the set of all vertices in the cloud
// close enough to 'ver', previously computed in 'relocate'
// we can destroy it here, it won't be used anymore
	
// see paragraph 9.8 in the manual

{	if ( not ver.exists() )  return false;  // no touch
	if ( not ver.belongs_to ( interface, tag::not_oriented ) )  return false;
	//  because 'ver' might have been left behind in the meanwhile

	assert ( set_of_ver.find(ver.core) == set_of_ver.end() );
	Cell prev_seg = Mesh::Progressive::interface.cell_behind(ver);
	Cell prev_ver = prev_seg.base().reverse();
	Cell next_seg = Mesh::Progressive::interface.cell_in_front_of(ver);
	Cell next_ver = next_seg.tip();
	set_of_ver.erase ( prev_ver.core );
	set_of_ver.erase ( next_ver.core );
	if ( set_of_ver.empty() )  return false;  // no touch
	if ( set_of_ver.size() == 1 )  return false;  // almost touch, no merge
	// do we need to check orientations ? probably not

	// we deal with two vertices or three vertices
	std::set<Cell::Core*>::iterator iit = set_of_ver.begin();
	assert ( iit != set_of_ver.end() );
	Cell one ( tag::whose_core_is, *iit );
	assert ( prev_ver != one );  assert ( next_ver != one );
	iit++;  assert ( iit != set_of_ver.end() );
	Cell two ( tag::whose_core_is, *iit );
	assert ( prev_ver != two );  assert ( next_ver != two );
	iit++;
	if ( iit != set_of_ver.end() )  // three vertices
	{	Cell three ( tag::whose_core_is, *iit );
		iit++;  assert ( iit == set_of_ver.end() );
		// we want one-two-three to be in order along the interface
		Cell seg = Mesh::Progressive::interface.cell_behind (one);
		Cell tmp = seg.base().reverse();
		if ( set_of_ver.find ( tmp.core ) != set_of_ver.end() )  // 'one' is not the first one
		{	if ( tmp == two )  { two = one;  one = tmp;  }
			else  {	assert ( tmp == three );  three = one;  one = tmp;  }             }
		seg = Mesh::Progressive::interface.cell_behind ( one );
		tmp = seg.base().reverse();
		if ( set_of_ver.find ( tmp.core ) != set_of_ver.end() )
		// 'one' is still not the first one
		{	if ( tmp == two )
			{ two = one;  one = tmp;  }
			else  {	assert ( tmp == three );  three = one;  one = tmp;  }  }
		// now 'one' must be the first one
		seg = Mesh::Progressive::interface.cell_behind ( one );
		assert ( set_of_ver.find ( seg.base().reverse().core ) == set_of_ver.end() );
		Cell one_two = Mesh::Progressive::interface.cell_in_front_of(one);
		two = one_two.tip();
		Cell two_three = Mesh::Progressive::interface.cell_in_front_of(two);
		three = two_three.tip();
		if ( Manifold::working.dist_sq ( one, next_ver ) < Mesh::Progressive::long_dist_sq )
		{	// std::cout << "touching interface S fill" << std::endl << std::flush;
			//  interface.cell_in_front_of(next_ver) may have tip 'one'
			Cell ver_two = Mesh::Progressive::glue_two_segs_S
				( ver, next_ver, next_seg, one, two, one_two );
			Cell ver_three ( tag::segment, ver.reverse(), three );
			Mesh::Progressive::fill_60 ( ver_two, two_three, ver_three.reverse(),
		                               two, Mesh::Progressive::interface         );
			ver_three.add_to ( Mesh::Progressive::interface );
			Mesh::Progressive::build_one_normal ( ver, three, ver_three );  // based on previous segment
			return true;                                                                }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < Mesh::Progressive::long_dist_sq )
		{	// std::cout << "touching interface Z fill" << std::endl << std::flush;
			//  interface.cell_behind(prev_ver) may have base 'three'
			Cell two_ver = Mesh::Progressive::glue_two_segs_Z
				( prev_ver, ver, prev_seg, two, three, two_three );
			Cell one_ver ( tag::segment, one.reverse(), ver );
			Mesh::Progressive::fill_60 ( one_two, two_ver, one_ver.reverse(),
		                               two, Mesh::Progressive::interface    );
			one_ver.add_to ( Mesh::Progressive::interface );
			Mesh::Progressive::build_one_normal ( one, ver, one_ver );  // based on previous segment
			return true;                                                                        }  }
	else  // two vertices
	{	// if 'two' is not next to 'one' within the interface, switch them
		Cell seg = Mesh::Progressive::interface.cell_in_front_of (one);
		if ( seg.tip() != two )
		{	// in some cases, 'one' and 'two' are not even adjacent to each other
			// sac20091971
			Cell other_seg = Mesh::Progressive::interface.cell_behind ( one );
			if ( other_seg.base().reverse() != two ) return false;  // no merge
			Cell tmp = one;  one = two;
			two = tmp;                            }
		Cell one_two = Mesh::Progressive::interface.cell_in_front_of(one);
		assert ( one_two.tip() == two );
		if ( Manifold::working.dist_sq ( one, next_ver ) < Mesh::Progressive::long_dist_sq )
		{	// std::cout << "touching interface S" << std::endl << std::flush;
			//  interface.cell_in_front_of(next_ver) may have tip 'one'
			Mesh::Progressive::glue_two_segs_S ( ver, next_ver, next_seg, one, two, one_two );
			return true;                                                                        }
		if ( Manifold::working.dist_sq ( two, prev_ver ) < Mesh::Progressive::long_dist_sq )
		{	// std::cout << "touching interface Z" << std::endl << std::flush;
			//  interface.cell_behind(prev_ver) may have base 'two'
			Mesh::Progressive::glue_two_segs_Z ( prev_ver, ver, prev_seg, one, two, one_two );
			return true;                                                                        }  }

	return false;  // almost touch, no merge
	
}  // end of Mesh::Progressive::check_touching

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::build_one_normal ( Cell & B, Cell & C, Cell & new_seg )  // static

// builds the normal vector for 'new_seg', based on information from previous segment
	
{	assert ( B == new_seg.base().reverse() );
	assert ( C == new_seg.tip() );
	Cell AB = Mesh::Progressive::interface.cell_behind ( B );
	assert ( AB.tip() == B );
	Cell A = AB.base().reverse();

	std::vector < double > vA = Manifold::working.coordinates() ( A );
	std::vector < double > vB = Manifold::working.coordinates() ( B );
	std::vector < double > old_e ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  old_e[i] = vB[i] - vA[i];
	assert ( Mesh::Progressive::normals.find(AB.core) != Mesh::Progressive::normals.end() );
	std::vector < double > & old_f = Mesh::Progressive::normals[AB.core];
	// 'e' is the vector of the segment, 'f' is orthogonal
	// they are all of approximately the same length, equal to Mesh::Progressive::desired_length

	// code below is identical to part of Mesh::Progressive::build_each_normal
	// if you change anything, please change both; keep them identical
	std::vector < double > vC = Manifold::working.coordinates() ( C );
	std::vector < double > new_e ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  new_e[i] = vC[i] - vB[i];
	// scalar products :
	double with_e = 0.,  with_f = 0., norm_e_sq = 0., norm_f_sq = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];     with_f += new_e[i]*old_f[i];
		norm_e_sq += old_e[i]*old_e[i];  norm_f_sq += old_f[i]*old_f[i];  }
	with_e /= norm_e_sq;
	with_f /= norm_f_sq;
	std::vector < double > new_f ( Mesh::Progressive::nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	Mesh::Progressive::improve_normal ( B, C, new_e, new_f );
	Mesh::Progressive::normals[new_seg.core] = new_f;  // optimize
	// code above is identical to part of Mesh::Progressive::build_each_normal
}  // end of Mesh::Progressive::build_one_normal

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::build_each_normal  // static
( Cell & B, Cell & C, Cell & new_seg, std::vector < double > & old_e, std::vector < double > & old_f )

// 'e' is the vector of the segment, 'f' is orthogonal
// they are all of approximately the same length, equal to Mesh::Progressive::desired_length
// the computed normal is 'new_f'
// at the end the old vectors 'e' and 'f' are replaced by the new ones

{	assert ( B == new_seg.base().reverse() );
	assert ( C == new_seg.tip() );
	std::vector < double > vB = Manifold::working.coordinates() ( B );
	// code below is identical to part of Mesh::Progressive::build_one_normal
	// if you change anything, please change both; keep them identical
	std::vector < double > vC = Manifold::working.coordinates() ( C );
	std::vector < double > new_e ( Mesh::Progressive::nb_of_coords );
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  new_e[i] = vC[i] - vB[i];
	// scalar products :
	double with_e = 0.,  with_f = 0.;
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )  // use Riemannian metric !
	{	with_e += new_e[i]*old_e[i];  with_f += new_e[i]*old_f[i];  }
	with_e /= Mesh::Progressive::desired_length_sq;
	with_f /= Mesh::Progressive::desired_length_sq;
	std::vector < double > new_f ( Mesh::Progressive::nb_of_coords );
	// we rotate 'new_e' with 90 degrees, in the same sense as 'f' is rotated from 'e'
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ )
		new_f[i] = - with_f * old_e[i] + with_e * old_f[i];
	// project and normalize :
	Mesh::Progressive::improve_normal ( B, C, new_e, new_f );
	Mesh::Progressive::normals[new_seg.core] = new_f;  // optimize
	// code above is identical to part of Mesh::Progressive::build_one_normal
	old_e = new_e;  old_f = new_f;
}  // end of Mesh::Progressive::build_each_normal

//-------------------------------------------------------------------------------------------------

inline Cell Mesh::Progressive::build_normals ( const Cell & start )  // static

// from a cell 'start', propagate normals along Mesh::Progressive::interface
// (will only cover the connected component containing 'start')
// return the first segment which already has a normal
// see paragraph 9.6 in the manual
// 'normal' should have norm approximately equal to Mesh::Progressive::desired_length

{	assert ( start.belongs_to ( Mesh::Progressive::interface, tag::oriented ) );
	Cell seg = start;
	// std::cout << "propagating normals" << std::endl << std::flush;
	Cell A = seg.base().reverse(),  B = seg.tip();
	std::vector < double > va = Manifold::working.coordinates() (A),
		vb = Manifold::working.coordinates() (B);
	assert ( Mesh::Progressive::normals.find ( seg.core ) != Mesh::Progressive::normals.end() );
	std::vector < double > e ( Mesh::Progressive::nb_of_coords ),
		f = Mesh::Progressive::normals[seg.core];  // assert key found
	for ( size_t i = 0; i < Mesh::Progressive::nb_of_coords; i++ ) e[i] = vb[i] - va[i];
	// 'e' and 'f' form an oriented basis in the two-dimensional space
	// tangent to the manifold at the current point
	// they will be used to build further vectors pointing outwards
	// (on the correct side of Mesh::Progressive::interface)
	while ( true )
	// Mesh::Progressive::interface may be disconnected, so we cannot use CellIterators
	// this loop will only cover its current connected component
	{ Cell new_seg = Mesh::Progressive::interface.cell_in_front_of ( B );
		if ( Mesh::Progressive::normals.find(new_seg.core) !=
				 Mesh::Progressive::normals.end()                   )  return new_seg;
		assert ( new_seg != start );
		Cell C = new_seg.tip();
		Mesh::Progressive::build_each_normal ( B, C, new_seg, e, f );
		// 'e' and 'f' get updated within 'build_each_normal'
		seg = new_seg;  B = C;                                                       }
	
}  // end of Mesh::Progressive::build_normals

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::glue_two_segs_common  // static
( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD, Cell & AD, Cell & BC )

{	AB.remove_from ( Mesh::Progressive::interface );
	CD.remove_from ( Mesh::Progressive::interface );
	Cell CB = BC.reverse();
	AD.add_to ( Mesh::Progressive::interface );
	CB.add_to ( Mesh::Progressive::interface );
	assert ( Mesh::Progressive::normals.find(AB.core) != Mesh::Progressive::normals.end() );
	Mesh::Progressive::normals.erase ( AB.core );  // assert that key was found
	assert ( Mesh::Progressive::normals.find ( CD.core ) != Mesh::Progressive::normals.end() );
	Mesh::Progressive::normals.erase ( CD.core );
	// build normals for two newly added segments AD and CB
	Mesh::Progressive::build_one_normal ( A, D, AD );  // from previous segment
	Mesh::Progressive::build_one_normal ( C, B, CB );  // from previous segment
}  // end of void Mesh::Progressive::glue_two_segs_common

//-------------------------------------------------------------------------------------------------

inline Cell Mesh::Progressive::glue_two_segs_S  // static
( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD )

// see paragraph 9.8 in the manual

// interface.cell_in_front_of(B) may have tip C
// that is, BC may belong already to 'interface'

{ Cell AD ( tag::segment, A.reverse(), D );
	Cell AC ( tag::segment, A.reverse(), C );
	Cell BC = interface.cell_in_front_of(B);
	if ( BC.tip() == C )
	{	Mesh::Progressive::fill_60 ( AB, BC, AC.reverse(), B, interface );
		AC.add_to ( interface );
		Mesh::Progressive::build_one_normal ( A, C, AC );  // based on previous segment
		Mesh::Progressive::fill_60 ( AC, CD, AD.reverse(), C, interface );
		AD.add_to ( interface );
		Mesh::Progressive::build_one_normal ( A, D, AD );  }  // based on previous segment
	else
	{	BC = Cell ( tag::segment, B.reverse(), C );
		Cell ABC ( tag::triangle, AB, BC, AC.reverse() );
		Cell ACD ( tag::triangle, AC, CD, AD.reverse() );
		ABC.add_to ( Mesh::Progressive::under_constr );
		ACD.add_to ( Mesh::Progressive::under_constr );
		Mesh::Progressive::glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return AD;                                                                    }

//-------------------------------------------------------------------------------------------------

inline Cell Mesh::Progressive::glue_two_segs_Z  // static
( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD )

// see paragraph 9.8 in the manual

// interface.cell_behind(A) may have base D
// that is, DA may belong already to 'interface'

{	Cell BC ( tag::segment, B.reverse(), C );
	Cell DB ( tag::segment, D.reverse(), B );
	Cell DA = interface.cell_behind(A);
	if ( DA.base().reverse() == D )
	{	Mesh::Progressive::fill_60 ( DA, AB, DB.reverse(), A, interface );
		DB.add_to ( interface );
		Mesh::Progressive::build_one_normal ( D, B, DB );  // based on previous segment
		Mesh::Progressive::fill_60 ( CD, DB, BC, D, interface );
		BC.reverse().add_to ( interface );
		Mesh::Progressive::build_one_normal ( B, C, BC );  }  // based on previous segment
	else
	{	DA = Cell ( tag::segment, D.reverse(), A );
		Cell ABD ( tag::triangle, AB, DB.reverse(), DA );
		Cell BCD ( tag::triangle, BC, CD, DB );
		ABD.add_to ( Mesh::Progressive::under_constr );
		BCD.add_to ( Mesh::Progressive::under_constr );
		Cell AD = DA.reverse();
		Mesh::Progressive::glue_two_segs_common ( A, B, AB, C, D, CD, AD, BC );  }
	return BC.reverse();                                                          }

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::fill_60  // static
( Cell & AB, Cell & BC, const Cell & CA, const Cell & B, Mesh & interface )

{	AB.remove_from ( interface );
	BC.remove_from ( interface );
	assert ( Mesh::Progressive::normals.find(AB.core) != Mesh::Progressive::normals.end() );
	assert ( Mesh::Progressive::normals.find(BC.core) != Mesh::Progressive::normals.end() );
	Mesh::Progressive::normals.erase ( AB.core );
	Mesh::Progressive::normals.erase ( BC.core );
	Mesh::Cloud::remove_point ( B );
	Cell new_tri ( tag::triangle, AB, BC, CA );
	new_tri.add_to ( under_constr );
	Mesh::Progressive::under_constr.baricenter ( B, AB );
}

//-------------------------------------------------------------------------------------------------

inline void Mesh::Progressive::fill_last_triangle ( const Cell & A, const Cell & B, const Cell & C,
	Cell & AB, Cell & BC, Cell & CA, Mesh & interface )  // static

{	Mesh::Progressive::fill_60 ( AB, BC, CA, B, interface );
	CA.remove_from ( interface );
	assert ( Mesh::Progressive::normals.find(CA.core) != Mesh::Progressive::normals.end() );
	Mesh::Progressive::normals.erase ( CA.core );
	Mesh::Cloud::remove_point ( A );
	Mesh::Cloud::remove_point ( C );
	Mesh::Progressive::under_constr.baricenter ( A, CA );
	Mesh::Progressive::under_constr.baricenter ( C, BC );                   }

//-------------------------------------------------------------------------------------------------

size_t Mesh::Progressive::nb_of_coords;

const double Mesh::Progressive::one_plus_tolerance = 1.2;
const double Mesh::Progressive::sqrt_of_075 = std::sqrt ( 0.75 );
const double Mesh::Progressive::half_of_sqrt_of_075 = Mesh::Progressive::sqrt_of_075 / 2.;
double Mesh::Progressive::desired_length, Mesh::Progressive::desired_length_sq,
	Mesh::Progressive::short_dist, Mesh::Progressive::short_dist_sq,
	Mesh::Progressive::long_dist, Mesh::Progressive::long_dist_sq;

Mesh Mesh::Progressive::under_constr ( tag::of_dimension_one );
Mesh Mesh::Progressive::interface ( tag::of_dimension_one );
// whatever; these variables will be set when we start the meshing process,
// in the body of the progressive Mesh constructor

//-------------------------------------------------------------------------------------------------

int main5 ( )

// skeleton of a tetrahedron
	
{	Manifold RR3 ( tag::Euclid, tag::of_dim, 3 );
	Function xyz = RR3.build_coordinate_system ( tag::Lagrange, tag::of_degree, 1 );
	Function x = xyz[0],  y = xyz[1],  z = xyz[2];

	// A ( 1, 0, 1 )
	// B ( 0, 1,-1 )
	// C (-1, 0, 1 )
	// D ( 0,-1,-1 )
	Function lambda_AB = ( x - y + 2.*z - 3. ) / 6.;
	Function dist_AB_sq = ( x - 1. - lambda_AB ) * ( x - 1. - lambda_AB )
		+ ( y + lambda_AB ) * ( y + lambda_AB )
		+ ( z - 1. - 2.*lambda_AB ) * ( z - 1. - 2.*lambda_AB );
	Function dist_AC_sq = y*y + (z-1.)*(z-1.);
	Function lambda_AD = ( x + y + 2.*z - 3. ) / 6.;
	Function dist_AD_sq = ( x - 1. - lambda_AD ) * ( x - 1. - lambda_AD )
		+ ( y - lambda_AD ) * ( y - lambda_AD )
		+ ( z - 1. - 2.*lambda_AD ) * ( z - 1. - 2.*lambda_AD );
	Function lambda_BC = ( - x - y + 2.*z + 3. ) / 6.;
	Function dist_BC_sq = ( x + lambda_BC ) * ( x + lambda_BC )
		+ ( y - 1. + lambda_BC ) * ( y - 1. + lambda_BC )
		+ ( z + 1. - 2.*lambda_BC ) * ( z + 1. - 2.*lambda_BC );
	Function dist_BD_sq = x*x + (z+1.)*(z+1.);
	Function lambda_CD = ( - x + y + 2.*z - 3. ) / 6.;
	Function dist_CD_sq = ( x + 1. + lambda_CD ) * ( x + 1. + lambda_CD )
		+ ( y - lambda_CD ) * ( y - lambda_CD )
		+ ( z - 1. - 2.*lambda_CD ) * ( z - 1. - 2.*lambda_CD );
	Function smd = smooth_min ( dist_AB_sq, dist_AC_sq, dist_AD_sq,
	                            dist_BC_sq, dist_BD_sq, dist_CD_sq, tag::threshold, 0.02 );
	Manifold tetra = RR3.implicit ( smd == 0.02 );
	// the product of distances could be a good example for non-uniform mesh generation
	// ( dist_AB_sq * dist_AC_sq * dist_AD_sq * dist_BC_sq * dist_BD_sq * dist_CD_sq == 0.2 );

	Cell A ( tag::vertex );  x(A) =  0.11;  y(A) = 0.   ;  z(A) = 1.2;
	Cell B ( tag::vertex );  x(B) =  0.09;  y(B) = 0.   ;  z(B) = 1.2;
	Cell C ( tag::vertex );  x(C) =  0.1 ;  y(C) = 0.017;  z(C) = 1.197;
	tetra.project(A);  tetra.project(B);  tetra.project(C);
	Cell AB ( tag::segment, A.reverse(), B );
	Cell BC ( tag::segment, B.reverse(), C );
	Cell CA ( tag::segment, C.reverse(), A );
	Mesh chain ( tag::of_dim_one );
	AB.add_to ( chain );
	BC.add_to ( chain );
	CA.add_to ( chain );

	std::vector < double > N { 0, -0.02, 0. };
	Mesh msh ( tag::progressive, tag::boundary, chain,
	           tag::start_at, AB, tag::towards, N,
	           tag::desired_length, 0.02              );

	msh.export_msh ("msh-new.msh");
	std::cout << "reached end, " << chain.core->cells[1].size() << " segments on the interface" << std::endl;

	return 0;
}

