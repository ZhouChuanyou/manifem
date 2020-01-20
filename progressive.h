
// manifem/progressive.h 2020.01.10

#ifndef MANIFEM_PROGRESSIVE_H
#define MANIFEM_PROGRESSIVE_H

#include <set>
#include <map>
#include <vector>
#include <forward_list>
#include <cmath>
#include <memory>
#include <iostream>
#include <sstream>
#include "assert.h"

#include "mesh.h"
#include "manifold.h"

namespace maniFEM {

	
// We have a cloud of points in a metric space, with immediate neighbours at distance approximately 1m.
// We choose "levels", like military ranks.
// To fix ideas, suppose we have "officers" of 10m, of 100m and of 1km.
// Direct subalterns of an officer of rank 'r' will be points closer than Cloud::distance[r]
// from that officer (not necessarily all those points).
// Indirect subalterns of an officer of rank 'r' will be points closer than Cloud::range[r]
// from that officer (not necessarily all of those).
// The union is the whole set. That is, for any given rank, all subordinates of all officers of that rank
// cover the entire set of points.
// Each vertex is represented by one or more officers.
// If it is only one, it has rank 0. Otherwise, ranks are 0, 1, 2, ..., rmax.

// Example
// We create one vertex only, it will become rank 0 officer.
// We add a second vertex. Both will be rank 0 officers, one of them will be also rank 1 officer.
// Actually, if the distance is large, it will also have versions of higher rank.
// We add a third vertex. Depends on the distances, but it will often be rank 0 officer,
// subaltern of existing rank 1 officer.

// Two officers of rank 'r' are always placed at distance > Cloud::distance[r]
// (otherwise, one of them should be subordinate of the other).

// Is this a quad-tree ?


namespace tag

{	struct BuildLowerRank { };  static const BuildLowerRank build_lower_rank;
	struct MinDist { };  static const MinDist min_dist;
	struct Ratio { };  static const Ratio ratio;
	struct OrthogonalTo { };  static const OrthogonalTo orthogonal_to;
	struct To { };  static const To to;                                        }


inline double approx_sqrt ( double arg, const tag::Around &, double centre, double sqrt_of_centre )

// we assume all segments have length around the same value, here called 'sqrt_of_centre'
// so the result should not be too far from 'sqrt_of_centre'
	
{	double res = 0.5 / sqrt_of_centre * ( arg + centre );
	#ifndef NDEBUG
	if ( ( res < 0.25*sqrt_of_centre ) or ( res > 4.*sqrt_of_centre ) )
	{	std::cout << "bad approximation of square root" << std::endl;
		std::cout << arg << " " << centre << std::endl;
		std::cout << res << " " << sqrt_of_centre << std::endl;
		exit(1);                                                         }
	#endif  // debug
	return res;
}


class Mesh::Officer

// a Point in the cloud, with a certain rank

{	public :

	class RankZero;  class PositiveRank;

	Cell point;

	// perhaps we should define 'rank' as a data member of Mesh::Officer
	// and eliminate 'get_rank', or define it inline
	
	Mesh::Officer::PositiveRank * superior { nullptr };
  // superior == nullptr if this is highest officer in Cloud

	inline Officer ( Cell p ) : point { p }
	{ assert ( p.dim() == 0 ); }

	virtual ~Officer ( )  { };

	inline Officer * promote ( size_t r );
	// promote 'this' to rank r, creating the necessary intermediates

	bool try_to_cling_to ( Mesh::Officer::PositiveRank * insert_place, Mesh::Officer * to_avoid );
	
	virtual size_t get_rank ( ) = 0;

	static std::map < Cell::Core *, Mesh::Officer::RankZero * > officer;  // use Cell::hook, it's faster

	virtual void find_close_vertices ( Cell cll, double dist_high, std::set < Cell::Core* > & res ) = 0;
	// among subalterns of 'this', returns all vertices closer to cll than dist_high

};  // end of class Mesh::Officer


class Mesh::Officer::RankZero : public Mesh::Officer

{	public :

	inline RankZero ( Cell p ) : Officer ( p )
	{	Mesh::Officer::officer[p.core] = this;  }
	
	size_t get_rank ( )  // virtual from Officer
	{	return 0;  }

	void find_close_vertices ( Cell cll, double dist_high, std::set < Cell::Core* > & res );
	// virtual from Mesh::Officer
	// among subalterns of 'this', returns all vertices closer to cll than dist_high
};


class Mesh::Officer::PositiveRank : public Mesh::Officer

{	public :

	size_t rank;
	// perhaps we should define 'rank' as a data member of Mesh::Officer
	// and eliminate 'get_rank', or define it inline

	std::set < Mesh::Officer * > subalterns;  // vertices in the neighbourhood

	inline PositiveRank ( Cell p, size_t r ) : Mesh::Officer ( p ), rank { r }  { }

	// receive a Cell (a vertex), create an officer on this point and
	// all lower rank copies of itself as the only subalterns
	inline PositiveRank ( Cell p, size_t r, const tag::BuildLowerRank & );

	size_t get_rank ( )  // virtual from Officer
	{	return rank;  }

	void adopt ( Mesh::Officer * of, Mesh::Officer * to_avoid = nullptr );

	inline void direct_adopt ( Mesh::Officer * of, size_t rank_of_this );
	
	void find_close_vertices ( Cell cll, double dist_high, std::set < Cell::Core* > & res );
	// virtual from Mesh::Officer
	// among subalterns of 'this', returns all vertices closer to cll than dist_high

	void print_subalterns ( size_t );
	
};  // end of class Mesh::Officer::PositiveRank


Mesh::Officer::PositiveRank::PositiveRank ( Cell p, size_t r, const tag::BuildLowerRank & )

// ('p' is a point) create an Officer on 'p' and
// all lower rank copies of itself as the only subalterns

:	Mesh::Officer::PositiveRank ( p, r )

{	assert ( r > 0 );
  Mesh::Officer * sub;
	if ( r == 1 ) sub = new Mesh::Officer::RankZero ( p );
	else sub = new Mesh::Officer::PositiveRank ( p, r-1, tag::build_lower_rank );
	std::pair < std::set<Mesh::Officer*>::iterator, bool > pa = this->subalterns.insert ( sub );
	assert ( pa.second );
	sub->superior = this;                                                                         }


class Mesh::Cloud

// essentially, Cloud is used as a namespace
// at most one Cloud object may exist
// no two meshing process may run simultaneously

{	public :

	// Manifold *manifold;

	static Mesh::Officer * highest_officer;

	static double ratio_between_distances;   // must be >2; 10 for example, or 7
	static std::vector < double > distance;  // size >= highest_rank + 1
	static std::vector < double > range;     // size == distance.size() >= highest_rank + 1
	// for example, distance could be [ 2, 14,  98, 686 ]
	// in that case, 'range' would be [ 0, 14, 112, 798 ]
	// 'distance' means :
	//    any direct subaltern of an officer is at distance <= distance[rank] of that officer
	// 'range' means :
	//    any indirect subaltern of an officer is at distance <= range[rank] of that officer
	// in the above, 'rank' is the rank of the higher officer, of course

	// constructors

	inline Cloud ( ) = delete;
	inline Cloud ( const Cloud & ) = delete;
	inline Cloud ( const Cloud && ) = delete;
	
	inline Cloud ( const tag::MinDist &, double d, const tag::Ratio &, double r )

	{	assert ( Mesh::Cloud::distance.size() == 0 );
		assert ( Mesh::Cloud::range.size() == 0 );
		// insert more explicit error message "at most one instance of Cloud may exist at a given time"
		Mesh::Cloud::distance.push_back ( d );
		Mesh::Cloud::range.push_back ( 0. );
		Mesh::Cloud::ratio_between_distances = r;
		Mesh::Cloud::highest_officer = nullptr;        }

	inline ~Cloud ( )
	{	Mesh::Cloud::distance.clear();
		Mesh::Cloud::range.clear();     }

	inline size_t highest_rank ()
	{	assert ( highest_officer );
		return highest_officer->get_rank();   }

	static void add_point ( Cell p );

	static void remove_point ( Cell p );
	
	static inline std::set<Cell::Core*> find_close_vertices ( Cell cll, double dist_high )
	// returns all vertices in the cloud closer to cll than dist_high
	{	std::set < Cell::Core * > res;
		Mesh::Cloud::highest_officer->find_close_vertices ( cll, dist_high, res );
		return res;                                                                  }
	
	static void print_everything ();

	static inline void build_vectors ( size_t r );
	// check and, if necessary, build vectors Cloud::distance and Cloud::range
};


inline Mesh::Officer * Mesh::Officer::promote ( size_t r )
// promote 'this' to rank r, creating the necessary intermediates
{	size_t this_rank = this->get_rank();
	assert ( r >= this_rank );
	assert ( Cloud::distance.size() > r );
	if ( r == this_rank ) return this;
	Mesh::Officer * of = this;
	for ( size_t i = this_rank+1; i <= r; i++ )
	{	of->superior = new Mesh::Officer::PositiveRank ( point, i );
		std::pair < std::set<Officer*>::iterator, bool > pa =
			of->superior->subalterns.insert ( of );
		assert ( pa.second );
		of = of->superior;                                        }
	return of;
}  // end of Officer::promote


inline void Mesh::Cloud::build_vectors ( size_t r )  // static
// check and, if necessary, build vectors Cloud::distance and Cloud::range
{	size_t old_size = Mesh::Cloud::distance.size();
	if ( old_size <= r )
	{	double d = Mesh::Cloud::distance.back();
		double last_range = Mesh::Cloud::range.back();
		Mesh::Cloud::distance.resize ( r+1 );
		Mesh::Cloud::range.resize ( r+1 );
		for ( size_t i = old_size; i <= r; i++ )
		{	d *= Cloud::ratio_between_distances;
			Mesh::Cloud::distance[i] = d ;
			last_range += d;
			Mesh::Cloud::range[i] = last_range;          }  }
}  // end of Cloud::build_vectors


inline void Mesh::Officer::PositiveRank::direct_adopt ( Mesh::Officer * of, size_t rank_of_this )

// we provide 'rank_of_this' just to avoid one more call to 'get_rank'

{	assert ( of );
	assert ( Mesh::Cloud::distance.size() > rank_of_this );
	assert ( Mesh::Cloud::distance.size() == Mesh::Cloud::range.size() );
	assert ( rank_of_this > of->get_rank() );
	// first promote 'of' then adopt it
	of = of->promote ( rank_of_this - 1 );
	of->superior = this;
	std::pair < std::set<Mesh::Officer*>::iterator, bool > pa = this->subalterns.insert ( of );
	assert ( pa.second );                                                                        }


struct Mesh::Progressive {  // used as a namespace

	static size_t nb_of_coords;
	
	static double desired_length, desired_length_sq;
	static const double one_plus_tolerance;
	static const double sqrt_of_075;
	static const double half_of_sqrt_of_075;
	static double short_dist, short_dist_sq, long_dist, long_dist_sq;

	// we are building a mesh 'under_constr'
	// we are using a moving interface which may be disconnected
	// (so we cannot use iterators on 'interface')
	static Mesh under_constr, interface;

	static std::map < Cell::Core *, std::vector < double > > normals;  // use Cell::hook

	// we use a temporary_vertex for projecting normals
	static Cell temporary_vertex;
	
	static inline bool positive_orientation
	( const Cell & A, const Cell & B, const Cell & AB, const Cell & CD );

	static inline double cos_sq_60
	( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & CD );

	static inline double cos_sq_120
	( const Cell & A, const Cell & B, const Cell & C, const Cell & AB, const Cell & CD );

	static inline void improve_normal
	( const Cell & A, const Cell & B, std::vector<double> & AB_coord, std::vector<double> & nor );

	static inline void build_one_normal ( Cell & B, Cell & C, Cell & new_seg );

	static inline void build_each_normal
	( Cell & B, Cell & C, Cell & new_seg, std::vector < double > & old_e,
	  std::vector < double > & old_f                                      );

	static inline Cell build_normals ( const Cell & start );

	static void advance ( const Cell & start );

	static inline void fill_60
		( Cell & AB, Cell & BC, const Cell & CA, const Cell & B, Mesh & interface );

	static inline void fill_last_triangle ( const Cell & A, const Cell & B, const Cell & C,
	  Cell & AB, Cell & BC, Cell & CA, Mesh & interface );

	static void relocate
	( const Cell & P, size_t n, std::vector<double> & direction,
	  const size_t & nc, std::set<Cell::Core*> & set_of_ver      );
	
	static inline bool check_touching ( Cell & P, std::set<Cell::Core*> & set_of_ver );
	
	static inline Cell glue_two_segs_S
	( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD );

	static inline Cell glue_two_segs_Z
	( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD );

	static inline void glue_two_segs_common
	( Cell & A, Cell & B, Cell & AB, Cell & C, Cell & D, Cell & CD, Cell & AD, Cell & BC );

};  // end of  struct Mesh::Progressive


}  // end of  namespace maniFEM

#endif
// ifndef MANIFEM_PROGRESSIVE_H
