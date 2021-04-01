
// mesh.h 2021.04.01

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


#ifndef MANIFEM_MESH_H
#define MANIFEM_MESH_H

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include "assert.h"

namespace maniFEM {

namespace tag {  // see paragraph 9.2 in the manual
	struct IsNegative { };  static const IsNegative is_negative;
	struct IsPositive { };  static const IsPositive is_positive;
	struct ThisMeshIsPositive { };  static const ThisMeshIsPositive this_mesh_is_positive;
	struct ReverseOrder { };  static const ReverseOrder reverse_order;
	struct ReverseOf { };  static const ReverseOf reverse_of;
	struct NonExistent { };  static const NonExistent non_existent;
	struct BuildIfNotExists { };  static const BuildIfNotExists build_if_not_exists;
	struct SeenFrom { };  static const SeenFrom seen_from;
	struct MayNotExist { };  static const MayNotExist may_not_exist;
	struct SurelyExists { };  static const SurelyExists surely_exists;
	struct CellsSurelyExist { };  static const CellsSurelyExist cells_surely_exist;
	struct OfDimension { };  static const OfDimension of_dim;
	                         static const OfDimension of_dimension;
	struct OfDimensionOne { };  static const OfDimensionOne of_dim_one;
	                            static const OfDimensionOne of_dimension_one;
	struct MinusOne { };  static const MinusOne minus_one;
	struct GreaterThanOne { };  static const GreaterThanOne greater_than_one;
	struct MightBeOne { };  static const MightBeOne might_be_one;
	struct Oriented { };  static const Oriented oriented;
	struct NotOriented { };  static const NotOriented not_oriented;
	struct DeepCopy { };  static const DeepCopy deep_copy;
	                      static const DeepCopy deep_copy_of;
	struct BuildCellsIfNec { };  static const BuildCellsIfNec build_cells_if_necessary;
	struct Progressive { };  static const Progressive progressive;
	struct StartAt { };  static const StartAt start_at;
	struct StopAt { };  static const StopAt stop_at;
	struct Towards { };  static const Towards towards;
	struct Boundary { };  static const Boundary boundary;
	struct BoundaryOf { };  static const BoundaryOf boundary_of;
	struct Bizarre { };  static const Bizarre bizarre;
	struct SizeMeshes { };  static const SizeMeshes size_meshes;
	struct BehindFace { };  static const BehindFace behind_face;
	struct InFrontOfFace { };  static const InFrontOfFace in_front_of_face;
	struct WithinMesh { };  static const WithinMesh within_mesh;
	struct OverCellsOfDim { };  static const CellsOfDim over_cells_of_dim;
	struct OverCellsOfReverseOf { }; static const CellsOfReverseOf over_cells_of_reverse_of;
	struct OverVertices { };  static const OverVertices over_vertices;
	struct OverSegments { };  static const OverSegments over_segments;
	struct Vertex { };  static const Vertex vertex; static const Vertex point;
	struct Segment { };  static const Segment segment;
	struct DividedIn { };  static const DividedIn divided_in;
	struct Triangle { };  static const Triangle triangle;
	struct Quadrangle { };  static const Quadrangle quadrangle;
	                        static const Quadrangle rectangle;
	                        static const Quadrangle square;
	                        static const Quadrangle quadrilateral;
	struct WhoseBoundaryIs { };  static const WhoseBoundaryIs whose_bdry_is;
	                             static const WhoseBoundaryIs whose_boundary_is;
	struct WhoseCoreIs { };  static const WhoseCoreIs whose_core_is;
	struct ForcePositive { };  static const ForcePositive force_positive;
	struct HasSize { };  static const HasSize has_size;
	struct ReserveSize { };  static const ReserveSize reserve_size;
	struct Pretty { };  static const Pretty pretty;
	struct Adapt { };  static const Adapt adapt;
	struct OfDegree { };  static const OfDegree of_degree;
	struct MeshIsBdry { };  static const MeshIsBdry mesh_is_bdry;
	struct MeshIsNotBdry { };  static const MeshIsNotBdry mesh_is_not_bdry;
	enum WithTriangles { with_triangles, not_with_triangles };
	struct Join { };  static const Join join;
	enum KeyForHook { tangent_vector, normal_vector, node_in_cloud };
	struct Onto { };  static const Onto onto;
	struct EntireManifold { };  static const EntireManifold entire_manifold;
	struct DesiredLength { };  static const DesiredLength desired_length;
	struct IntrinsicOrientation { };  static const IntrinsicOrientation intrinsic_orientation;
	struct InherentOrientation { };  static const InherentOrientation inherent_orientation;
	struct RandomOrientation { };  static const RandomOrientation random_orientation;             }

class Cell;  class Mesh;
class CellIterator;  class MeshIterator;
class Manifold;  class Function;

//-----------------------------------------------------------------------------//


namespace Util {

inline size_t assert_diff ( size_t a, size_t b )
{	assert ( a >= b );  return  a - b;  }

template < typename X, typename Y > inline Y assert_cast ( X x )
#ifndef NDEBUG
{	Y y = dynamic_cast < Y > (x);  assert (x);  }
#else
{	Y y = static_cast < Y > (x);  }
#endif

}  // namespace maniFEM::Util


//-----------------------------------------------------------------------------//
//-----------------  wrappers Cell and Mesh  ----------------------------------//
//-----------------------------------------------------------------------------//


// a cell of dimension zero is a point, see class Cell::Positive::Vertex and Negative::Vertex
// a cell of dimension one is a segment, see class Cell::Positive::Segment and Negative::Segment
// for cells of dimension two or more, see class Cell::Positive::HighDim and Cell::Negative::HighDim
// a cell of dimension two may be a triangle, a quadrangle or some other polygon
// a cell of dimension three may be a tetrahedron, a cube or some other polyhedron
// cells of dimension four or higher may be constructed,
// but their usefulness is questionable

// cells may be positively or negatively oriented
// see class Cell::Positive and class Cell::Negative
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// a cell is mainly defined by its boundary (which is a mesh of lower dimension)
// the orientation of a cell is nothing more than an orientation of its boundary
// see the comments on orientation in class Mesh below

// see paragraphs 1.2, 8.1, 8.7 and 9.3 in the manual

class Cell

// a thin wrapper around a Cell::Core, with most methods delegated to 'core'

{	public :

	class Core;
	
	Cell::Core * core;

	// many constructors defined after class Mesh::Negative

	inline Cell ( const tag::NonExistent & )
	:	core { nullptr }
	{ }

	inline Cell ( const tag::WhoseCoreIs &, Cell::Core * c )
	:	core { c }
	{	assert ( c );  }

	inline Cell ( const tag::WhoseBoundaryIs &, Mesh & );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell,
                const tag::BuildIfNotExists & build = tag::build_if_not_exists );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::SurelyExists & );
	inline Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::MayNotExist & );
	// inline Cell ( const tag::WhoseBoundaryIs &, Mesh::Core * );
	inline Cell ( const tag::Vertex &, const tag::IsPositive & ispos = tag::is_positive );
	inline Cell ( const tag::Segment &, const Cell & A, const Cell & B );
	inline Cell ( const tag::Bizarre &, const Cell &, const Cell & );
	inline Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA );
	inline Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                         const Cell & CD, const Cell & DA );
	inline Cell ( const tag::BehindFace &, const Cell &, const tag::WithinMesh &,
	              const Mesh &, const tag::SurelyExists & se = tag::surely_exists );
	inline Cell ( const tag::BehindFace &, const Cell &,
                const tag::WithinMesh &, const Mesh &, const tag::MayNotExist & );
	inline Cell ( const tag::InFrontOfFace &, const Cell &, const tag::WithinMesh &,
                const Mesh &, const tag::SurelyExists & se = tag::surely_exists     );
	inline Cell ( const tag::InFrontOfFace &, const Cell &,
                const tag::WithinMesh &, const Mesh &, const tag::MayNotExist & );
	// use : Cell tri ( tag::in_front_of_face, f, tag::within_mesh, msh );

	// we are still in class Cell
	
	// methods delegated to 'core'

	inline Cell reverse
	( const tag::BuildIfNotExists & build = tag::build_if_not_exists ) const;
	inline Cell reverse ( const tag::SurelyExists & ) const;
	inline Cell reverse ( const tag::MayNotExist & ) const;
	inline Mesh boundary () const;
	inline bool exists () const  { return core != nullptr;  }
	inline bool is_positive () const;
	inline Cell get_positive ();
	inline size_t dim () const;
	inline bool has_reverse () const;
	inline Cell tip () const;
	inline Cell base () const;

	inline bool belongs_to ( const Mesh & msh, const tag::Oriented & ) const;
	inline bool belongs_to ( const Mesh & msh, const tag::NotOriented & ) const;

	// method 'glue_on_bdry_of' is intensively used when building a mesh
	// it glues 'this' cell to the boundary of 'cll'
	inline void glue_on_bdry_of ( Cell & cll );
	
	// method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll' - used mainly in remeshing
	inline void cut_from_bdry_of ( Cell & cll );
	
	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// if 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead
	inline void add_to ( Mesh & msh );
	inline void remove_from ( Mesh & msh );

	// we are still in class Cell

	inline void project ( ) const;
	
	inline void project ( const tag::Onto &, const Manifold m ) const;
		
#ifndef NDEBUG
	inline void print_everything ( );
#endif

	void print_coords();

	static std::vector < size_t > double_heap_size_pos, double_heap_size_neg,
		size_t_heap_size_pos, size_t_heap_size_neg, short_int_heap_size_pos, short_int_heap_size_neg;

	struct field_to_meshes
	{	short int counter_pos;
		short int counter_neg;
		std::list<Cell::Core*>::iterator where;
		inline field_to_meshes ( short int i, short int j )
		:	counter_pos {i}, counter_neg {j} { }
		inline field_to_meshes ( short int i, short int j,
		                         std::list<Cell::Core*>::iterator w )
		:	counter_pos {i}, counter_neg {j}, where {w} { }              };

	struct field_to_meshes_same_dim
	{	short int sign;
		std::list<Cell::Core*>::iterator where;
		inline field_to_meshes_same_dim ( short int i ) : sign {i} { }
		inline field_to_meshes_same_dim
		( short int i, std::list<Cell::Core*>::iterator w )
		:	sign {i}, where {w} { }                                       };

  class Positive;  class Negative;
	
}; // end of  class Cell


inline bool operator== ( const Cell & c1, const Cell & c2 )
{	return c1.core == c2.core;  }

inline bool operator!= ( const Cell & c1, const Cell & c2 )
{	return c1.core != c2.core;  }

inline bool operator< ( const Cell & c1, const Cell & c2 )
{	return c1.core < c2.core;  }

//-----------------------------------------------------------------------------//


// roughly speaking, a mesh is a collection of cells of the same dimension

// an orientation of the mesh
// is nothing more than an orientation of each of its cells (of maximum dimension)
// but these orientations cannot be arbitrary, they must be compatible
// in the sense that a face common to two cells must be seen as positive
// from one of the cells and as negative from the other cell

// as a consequence, for connected meshes there are only two possible orientations
// although nothing prevents a mesh to be disconnected, we only allow for two orientations
// which is which depends only on the construction process
// the first one to be created is positive, the other one will be negative,
// created when we call the 'reverse' method

// negative meshes will appear mostly as boundaries of negative cells
// that's why we do not store their core in the computer's memory
// wrappers for negative meshes are built on-the-fly, e.g. in method Cell::boundary
// their core points to the (positive) reverse, the only difference is in 'sign'

// see paragraphs 1.2, 8.1, 8.7 and 9.3 in the manual

class Mesh

// a thin wrapper around a Mesh::Core, with most methods delegated to 'core'

{	public :
	
	class Core;

	Mesh::Core * core;  // there are no cores for negative meshes
	// instead, we keep here a pointer to the direct (positive) core
	// thus, for meshes, 'core' points always to a positive Mesh::Core

	class Sign;  // emulate virtual methods
	// this is how we distinguish between a positive mesh and a negative one
	const Sign sign;

	// we keep here the topological dimension of the largest mesh we intend to build
	// see method 'set_max_dim' and paragraph 9.5 in the manual
	static size_t maximum_dimension_plus_one;

	// constructors :
	
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *,
                const tag::IsPositive & ispos = tag::is_positive );
	// builds a negative mesh from a positive one, assuming all cells have reverse :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::IsNegative &,
	              const tag::CellsSurelyExist &                                    );
	// builds a negative mesh from a positive one, creating reverse cells if necessary :
	inline Mesh ( const tag::WhoseCoreIs &, Mesh::Core *, const tag::IsNegative &,
                const tag::BuildCellsIfNec &                                     );

	Mesh ( const tag::DeepCopy &, const Mesh & );

	// we are still in class Mesh
	// geometric constructors are defined in global.cpp
	// segment, triangle, quadrangle, join
	
	inline Mesh ( const tag::Segment &, const Cell & A, const Cell & B, const tag::DividedIn &, size_t n );
	inline Mesh ( const tag::Pretty &, const tag::Segment &,
                const Cell & A, const Cell & B, const tag::DividedIn &, size_t n );
	inline Mesh ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	inline Mesh ( const tag::Pretty &, const tag::Triangle &,
                const Mesh & AB, const Mesh & BC, const Mesh & CA );
	inline Mesh ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	                                       const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles         );
	inline Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	                                       const Cell & NE, const Cell & NW, size_t m, size_t n,
	              const tag::WithTriangles & wt = tag::not_with_triangles                        );
	inline Mesh ( const tag::Pretty &, const tag::Quadrangle &,
	              const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	              const tag::WithTriangles & wt = tag::not_with_triangles );
	inline Mesh ( const tag::Join &, const std::list<Mesh> & l );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );
	inline Mesh ( const tag::Join &,
	              const Mesh &, const Mesh &, const Mesh &, const Mesh &, const Mesh & );
 
	// we are still in class Mesh
	// constructors with tag::Progressive are defined in progressive.cpp
	
	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
         const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length,
	       const tag::RandomOrientation &                                                 );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
	       const tag::DesiredLength &, const Function & length, const tag::RandomOrientation & );

	Mesh ( const tag::Progressive &, const tag::DesiredLength &, const Function & length,
	       const tag::InherentOrientation &                                               );

	Mesh ( const tag::Progressive &, const tag::EntireManifold, Manifold manif,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length              );

	inline Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	              const tag::DesiredLength &, const Function & length,
	              const tag::RandomOrientation &                                   )
	:	Mesh ( tag::progressive, tag::boundary, interface, tag::desired_length, length )  { }

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length, const tag::IntrinsicOrientation & );
 
	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	Mesh ( const tag::Progressive &, const tag::Boundary &, Mesh interface,
	       const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> normal,
	       const tag::DesiredLength &, const Function & length               );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> tangent,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                 );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::Towards &, std::vector<double> tangent,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length                  );

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::DesiredLength &, const Function & length, const tag::RandomOrientation & );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::StopAt &, const Cell & stop,
	              const tag::DesiredLength &, const Function & length,
	              const tag::RandomOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, stop,
		       tag::desired_length, length                                 )  { }

	Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	       const tag::StopAt &, const Cell & stop,
	       const tag::DesiredLength &, const Function & length, const tag::InherentOrientation & );

	inline Mesh ( const tag::Progressive &, const tag::StartAt &, const Cell & start,
	              const tag::DesiredLength &, const Function & length,
	              const tag::InherentOrientation & )
	:	Mesh ( tag::progressive, tag::start_at, start, tag::stop_at, start,
		       tag::desired_length, length, tag::inherent_orientation       )  { }

	// we are still in class Mesh
	
	void pretty_constructor ( const tag::Segment &, const Cell & A, const Cell & B,
	                          const tag::DividedIn &, size_t n );
	void pretty_constructor ( const tag::Triangle &, const Mesh & AB,
	                          const Mesh & BC, const Mesh & CA         );
	void pretty_constructor ( const tag::Quadrangle &, const Mesh & south, const Mesh & east,
	             const Mesh & north, const Mesh & west, const tag::WithTriangles & wt        );

	inline bool is_positive () const;
	inline size_t dim () const;
	inline Mesh reverse () const;

	inline size_t number_of ( const tag::CellsOfDim &, size_t d ) const;
	inline size_t number_of ( const tag::Vertices & ) const;
	inline size_t number_of ( const tag::Segments & ) const;

	inline Cell first_vertex ( ) const;
	inline Cell last_vertex ( ) const;
	inline Cell first_segment ( ) const;
	inline Cell last_segment ( ) const;

	// we are still in class Mesh
	
	inline Cell cell_in_front_of
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell cell_behind ( const Cell face, const tag::MayNotExist & ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::MayNotExist & ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SurelyExists & se = tag::surely_exists ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::MayNotExist & ) const;
		
	// we are still in class Mesh
	// the eight methods below are only relevant for STSI meshes

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::SurelyExists & se = tag::surely_exists              ) const;

	inline Cell cell_in_front_of
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline Cell cell_behind
	( const Cell face, const tag::SeenFrom &, const Cell neighbour,
	  const tag::MayNotExist &                                       ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;

	inline Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::MayNotExist &                                                       ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;

	inline Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::MayNotExist &                                                       ) const;

	// we are still in class Mesh
	
	void join_list ( const std::list<Mesh> & l );

	// method baricenter defined in manifold.h
	void inline baricenter ( const Cell & ver, const Cell & seg );

	inline CellIterator iterator ( const tag::OverVertices & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator ( const tag::OverVertices &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator ( const tag::OverSegments & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::ForcePositive & ) const;
	inline CellIterator iterator ( const tag::OverSegments &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ReverseOrder & ) const;
	inline CellIterator iterator ( const tag::OverCellsOfDim &, size_t d ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::Reverse & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive & ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &                                        ) const;
	inline CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	  const tag::ReverseOrder &                                          ) const;

	// methods draw_ps and export_msh defined in global.cpp

	void draw_ps ( std::string file_name );
	void draw_ps_3d ( std::string file_name );
	void export_msh ( std::string f, std::map<Cell::Core*,size_t> & ver_numbering );
	void export_msh ( std::string f );
	
#ifndef NDEBUG
	inline void print_everything ( );
#endif
	
	inline static void set_max_dim ( size_t d )
	// see paragraph 9.5 in the manual
	{	maximum_dimension_plus_one = d + 1;
		Cell::double_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::double_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );
		Cell::size_t_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::size_t_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_pos.resize ( maximum_dimension_plus_one, 0 );
		Cell::short_int_heap_size_neg.resize ( maximum_dimension_plus_one, 0 );  }
	
	inline static void add_link
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	inline static void remove_link
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	inline static void add_link_rev
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );
	inline static void remove_link_rev
		( Cell::Core *, Cell::Core *, Mesh::Core *, short int, short int );

	struct Connected  {  class OneDim;  class HighDim;  };
	struct MultiplyConnected  {  class OneDim;  class HighDim; };
	class ZeroDim;  class Fuzzy;  class STSI;
	
	// here is where the low-level linking between cells and meshes happens
	// ***_deep_connections and ***_deep_connections_rev
	// are called from Mesh::***::add_*** and remove_***
	static inline void make_deep_connections
	( Cell::Positive::NotVertex * cll, Mesh::Core * msh, const tag::MeshIsNotBdry & );
	static inline void break_deep_connections
	( Cell::Positive::NotVertex * cll, Mesh::Core * msh, const tag::MeshIsNotBdry & );
	static inline void make_deep_connections_rev
  ( Cell::Core * o_cll, Cell::Positive::NotVertex * cll,
	  Mesh::Core * msh, const tag::MeshIsNotBdry &         );
	static inline void break_deep_connections_rev
	( Cell::Core * o_cll, Cell::Positive::NotVertex * cll,
	  Mesh::Core * msh, const tag::MeshIsNotBdry &         );
	static inline void make_deep_connections
	( Cell::Positive::NotVertex * cll, Mesh::Core * msh, const tag::MeshIsBdry & );
	static inline void break_deep_connections
	( Cell::Positive::NotVertex * cll, Mesh::Core * msh, const tag::MeshIsBdry & );
	static inline void make_deep_connections_rev
  ( Cell::Core * o_cll, Cell::Positive::NotVertex * cll,
	  Mesh::Core * msh, const tag::MeshIsBdry &            );
	static inline void break_deep_connections_rev
	( Cell::Core * o_cll, Cell::Positive::NotVertex * cll,
	  Mesh::Core * msh, const tag::MeshIsBdry &            );

}; // end of  class Mesh


inline bool operator== ( const Mesh & m1, const Mesh & m2 )
{	return m1.core == m2.core;  }


//-----------------------------------------------------------------------------//
//----------------      core Cells     ----------------------------------------//
//-----------------------------------------------------------------------------//


class Cell::Core

// abstract class
// specialized in Cell::Positive::{Vertex,Segment,HighDim} and
//                Cell::Negative::{Vertex,Segment,HighDim}

{	public :

	Cell::Core * reverse_p { nullptr };

	// we keep numeric values here :
	std::vector < double > double_heap;
	std::vector < size_t > size_t_heap;
	std::vector < short int > short_int_heap;
	// and heterogeneous information here :
	std::map < tag::KeyForHook, void * > hook;

	// if 'this' is a face of another cell and that other cell belongs to some mesh msh,
	// cell_behind_within[msh] keeps that cell
	// see methods Mesh::cell_behind and Mesh::cell_in_front_of
	std::map < Mesh::Core*, Cell::Core* > cell_behind_within;

	inline Core ( const tag::IsPositive &, const tag::OfDimension &, size_t d )
	:	double_heap ( Cell::double_heap_size_pos [d] ),
		size_t_heap ( Cell::size_t_heap_size_pos [d] ),
		short_int_heap ( Cell::short_int_heap_size_pos [d] )
	{ }

	inline Core ( const tag::IsNegative& ,
		const tag::ReverseOf &, Cell::Core * rev, const tag::OfDimension &, size_t d )
	:	reverse_p { rev },
		double_heap ( Cell::double_heap_size_neg [d] ),
		size_t_heap ( Cell::size_t_heap_size_neg [d] ),
		short_int_heap ( Cell::short_int_heap_size_neg [d] )
	{	assert ( rev );  }

	inline Core ( const tag::Ghost & )  { }

	virtual ~Core ( ) { };

	virtual bool is_positive ( ) const = 0;
	virtual Cell::Positive * get_positive ( ) = 0;
	virtual size_t get_dim ( ) const = 0;
	virtual Cell::Core * reverse ( const tag::BuildIfNotExists & ) = 0;
	virtual Cell::Core * build_reverse ( ) = 0;

	virtual Cell::Core * tip ();
	virtual Cell::Core * base ();
	// execution forbidden
	// overridden by Cell::Positive::Segment and Cell::Negative::Segment

	virtual Mesh boundary ( ) = 0;

	virtual bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const = 0;
	virtual bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const = 0;

	// Method 'glue_on_bdry_of' is intensively used when building a mesh,
	// e.g. within factory functions in Cell class.
	// It glues 'this' cell to the boundary of 'cll'.
	inline void glue_on_bdry_of ( Cell::Core * cll );

	// Method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll'. Used mainly in remeshing.
	inline void cut_from_bdry_of ( Cell::Core * cll );
	
	// the two methods below are only relevant for vertices
	// so we forbid execution for now and then override them in Cell::***tive::Vertex
  virtual void add_to_seg ( Cell::Positive::Segment * seg );
	virtual void remove_from_seg ( Cell::Positive::Segment * seg );

	// the four methods below are not relevant for vertices
  virtual void add_to_mesh ( Mesh::Core * msh ) = 0;
	virtual void remove_from_mesh ( Mesh::Core * msh ) = 0;
  virtual void add_to_bdry ( Mesh::Core * msh ) = 0;
	virtual void remove_from_bdry ( Mesh::Core * msh ) = 0;

	virtual void glue_on_my_bdry ( Cell::Core * ) = 0;
	virtual void cut_from_my_bdry ( Cell::Core * ) = 0;
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see above)

	virtual void dispose () { };
	// assert there are no meshes above
	// for each face, cut_from_bdry
	// if the face has no other meshes above, dispose of it
	// dispose the (now empty) boundary, then delete this cell

#ifndef NDEBUG
	std::string name;
	virtual std::string get_name ( ) = 0;
	virtual void print_everything ( ) = 0;
#endif

}; // end of class Cell::Core

//-----------------------------------------------------------------------------//


class Cell::Positive : public Cell::Core

// abstract class, introduces attribute  meshes
// and defines virtual methods  is_positive, get_positive, belongs_to, glue_common, cut_common

// specialized in Cell::Positive::{Vertex,Segment,HighDim}

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// the 'meshes' attribute keeps information about all meshes
	// "above" 'this' cell, that is, all meshes containing 'this' as a cell
	// it is indexed over the dimension of the mesh minus the dimension of 'this' 
	// for each mesh, it keeps a 'Cell::field_to_meshes' value, containing two counters
	// and an iterator into the 'cells' field of that particular mesh
	// the iterator is only meaningful for fuzzy or STSI meshes
	// of course this implies quite some amount of redundant information
	// but this redundancy makes the classes fast, especially for remeshing
	
	// since indices of vectors begin at zero
	// the keys of 'meshes[i]' will be meshes of dimension 'i + this->dim'
	// however, for meshes of the same dimension as 'this' cell,
	// we keep a different record in Cell::Positive::{Vertex,NotVertex}::meshes_same_dim
	// so meshes[0] will always be an empty map

	std::vector < std::map < Mesh::Core *, Cell::field_to_meshes > > meshes;

	inline Positive ( const tag::OfDimension &, size_t d,
                            const tag::SizeMeshes &, size_t sz  )
	:	Cell::Core ( tag::is_positive, tag::of_dim, d ),
		meshes ( sz )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell::Positive * get_positive ( );  // virtual from Cell::Core

	Cell::Core * reverse ( const tag::BuildIfNotExists & );  // virtual from Cell::Core

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	inline void glue_common ( Cell::Core * face );
	inline void cut_common ( Cell::Core * face );
	// do not use directly; called from glue_on_my_bdry and cut_from_my_bdry

#ifndef NDEBUG
	std::string get_name();
#endif
	
	class Vertex;  class NotVertex;  class Segment;  class HighDim;
	
}; // end of class Cell::Positive


class Cell::Negative : public Cell::Core

// abstract class, defines virtual methods
// is_positive, get_positive, reverse, build_reverse, belongs_to,
// glue_on_my_bdry, cut_from_my_bdry

// specialized in Cell::Negative::{Vertex,Segment,HighDim}

{	public :

	inline Negative ( const tag::OfDimension &, size_t d,
														const tag::ReverseOf &, Cell::Positive * rev )
	: Cell::Core ( tag::is_negative, tag::reverse_of, rev, tag::of_dim, d )
	{	}

	inline Negative ( const tag::Ghost & )
	: Cell::Core ( tag::ghost )
	{	}

	bool is_positive ( ) const;  // virtual from Cell::Core
	Cell::Positive * get_positive ( );  // virtual from Cell::Core
	// virtual size_t get_dim ( ) const = 0;  // declared in Cell::Core

	Cell::Core * reverse ( const tag::BuildIfNotExists & );  // virtual from Cell::Core
	Cell::Core * build_reverse ( );  // virtual from Cell::Core, execution forbidden

	bool belongs_to ( Mesh::Core *, const tag::Oriented & ) const;  // virtual from Cell::Core
	bool belongs_to ( Mesh::Core *, const tag::NotOriented & ) const;  // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

#ifndef NDEBUG
	std::string get_name();
#endif
	
	class Vertex;  class Segment;  class HighDim;
	
}; // end of class Cell::Negative

//-----------------------------------------------------------------------------//


class Cell::Positive::Vertex : public Cell::Positive

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// attribute  meshes  inherited from Cell::Positive

	// in 'meshes_same_dim' we keep record of meshes of dimension zero "above" 'this' vertex,
	// that is, of segments having 'this' extremity
	// the 'short int' is a sign, 1 or -1

	std::map < Cell::Positive::Segment *, short int > meshes_same_dim;
	
	inline Vertex ( );

	Vertex ( const Cell::Positive::Vertex & ) = delete;
	Vertex ( const Cell::Positive::Vertex && ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex & ) = delete;
	Cell::Positive::Vertex & operator= ( const Cell::Positive::Vertex && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	// reverse  defined by Cell::Positive
	Cell::Core * build_reverse ( );  // virtual from  Cell::Core
	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	// two versions of  belongs_to  defined by Cell::Positive

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// the two methods below are defined by Cell::Core but overriden here
	void add_to_seg ( Cell::Positive::Segment * seg ) override;
	void remove_from_seg ( Cell::Positive::Segment * seg ) override;

	// the six methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh );
	void glue_on_my_bdry ( Cell::Core * );
	void cut_from_my_bdry ( Cell::Core * );

	// glue_common  and  cut_common  defined by Cell::Positive
	
#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Positive::Vertex


class Cell::Negative::Vertex : public Cell::Negative

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	inline Vertex ( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p );

	Vertex ( const Cell::Negative::Vertex & ) = delete;
	Vertex ( const Cell::Negative::Vertex && ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex & ) = delete;
	Cell::Negative::Vertex & operator= ( const Cell::Negative::Vertex && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// reverse  defined by Cell::Negative
	// build_reverse  defined by Cell::Negative, execution forbidden
	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core, here execution forbidden

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core but overriden here
	void add_to_seg ( Cell::Positive::Segment * seg ) override;
	void remove_from_seg ( Cell::Positive::Segment * seg ) override;

	// the four methods below are virtual from Cell::Core, here execution forbidden
	void add_to_mesh ( Mesh::Core * msh );
	void remove_from_mesh ( Mesh::Core * msh );
  void add_to_bdry ( Mesh::Core * msh );
	void remove_from_bdry ( Mesh::Core * msh );

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Negative::Vertex

//-----------------------------------------------------------------------------//


class Cell::Positive::NotVertex : public Cell::Positive

// abstract class, useful only for introducing the attribute  meshes_same_dim
// specialized in Cell::Positive::{Segment,HighDim}

{	public :

	// in 'meshes_same_dim' we keep record of meshes "above" 'this' cell,
	//   of the same dimension
	// in Cell::field_to_meshes_same_dim, the 'short int sign' is a sign, 1 or -1
	// the iterator 'where' is only meaningful for fuzzy or STSI meshes

	std::map < Mesh::Core *, Cell::field_to_meshes_same_dim > meshes_same_dim;
	
}; // end of  class Cell::Positive::NotVertex

//-----------------------------------------------------------------------------//


class Cell::Positive::Segment : public Cell::Positive::NotVertex

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// attribute  meshes  inherited from Cell::Positive
	// attribute  meshes_same_dim  inherited from Cell::Positive::NotVertex

	Cell::Negative::Vertex * base_p { nullptr };
	Cell::Positive::Vertex * tip_p { nullptr };

	inline Segment ( Cell::Negative::Vertex * a, Cell::Positive::Vertex * b );

	Segment ( const Cell::Positive::Segment & ) = delete;
	Segment ( const Cell::Positive::Segment && ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment & ) = delete;
	Cell::Positive::Segment & operator= ( const Cell::Positive::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual from Cell::Core

	// reverse  defined by Cell::Positive
	Cell::Core * build_reverse ( );  // virtual from  Cell::Positive

	Cell::Core * tip () override; // virtual, overrides definition by Cell::Core
	Cell::Core * base () override; // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Positive::Segment


class Cell::Negative::Segment : public Cell::Negative

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	inline Segment ( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p );

	Segment ( const Cell::Negative::Segment & ) = delete;
	Segment ( const Cell::Negative::Segment && ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment & ) = delete;
	Cell::Negative::Segment & operator= ( const Cell::Negative::Segment && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core
	
	// reverse  defined by Cell::Negative
	// build_reverse  defined by Cell::Negative, execution forbidden
	Cell::Core * tip () override;  // virtual, overrides definition by Cell::Core
	Cell::Core * base () override;  // virtual, overrides definition by Cell::Core

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif
	
}; // end of  class Cell::Negative::Segment

//-----------------------------------------------------------------------------//


class Cell::Positive::HighDim : public Cell::Positive::NotVertex

// a cell of dimension >= 2

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	// attribute  meshes  inherited from Cell::Positive
	// attribute  meshes_same_dim  inherited from Cell::Positive::NotVertex

	Mesh::Core * boundary_p;

	inline HighDim ( const tag::OfDimension &, size_t d,
	                 const tag::WhoseBoundaryIs &, Mesh::Core * msh );
	inline HighDim ( const tag::WhoseBoundaryIs &, Mesh::Core * msh );
	inline HighDim ( const tag::Triangle &, Cell::Core * AB,
	                       Cell::Core * BC, Cell::Core * CA );
	inline HighDim ( const tag::Quadrangle &, Cell::Core * AB, Cell::Core * BC,
	                                          Cell::Core * CD, Cell::Core * DA );

	HighDim ( const Cell::Positive::HighDim & ) = delete;
	HighDim ( const Cell::Positive::HighDim && ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim & ) = delete;
	Cell::Positive::HighDim & operator= ( const Cell::Positive::HighDim && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Positive
	size_t get_dim ( ) const; // virtual

	// reverse  defined by Cell::Positive
	Cell::Core * build_reverse ( );  // virtual from  Cell::Core
	// tip  and  base  defined by Cell::Core, execution forbidden

	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Positive

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	void glue_on_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	void cut_from_my_bdry ( Cell::Core * ); // virtual from Cell::Core
	// we feel that 'glue_on_bdry_of' and 'cut_from_bdry_of' are more readable
	// so we suggest to use those (see class Cell::Core)

	// glue_common  and  cut_common  defined by Cell::Positive
	
#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

}; // end of  class Cell::Positive::HighDim


class Cell::Negative::HighDim : public Cell::Negative

// a negative cell of dimension >= 2

{	public :

	// attributes reverse_p, double_heap, size_t_heap, short_int_heap, hook
	// and cell_behind_within inherited from Cell::Core

	inline HighDim ( const tag::OfDimension, size_t d,
                        const tag::ReverseOf &, Cell::Positive * direct_cell_p );
	inline HighDim ( const tag::ReverseOf &, Cell::Positive * direct_cell_p );

	HighDim ( const Cell::Negative::HighDim & ) = delete;
	HighDim ( const Cell::Negative::HighDim && ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim & ) = delete;
	Cell::Negative::HighDim & operator= ( const Cell::Negative::HighDim && ) = delete;

	// is_positive  and  get_positive  defined by Cell::Negative
	size_t get_dim ( ) const; // virtual from Cell::Core

	// reverse  defined by Cell::Negative
	// build_reverse  defined by Cell::Negative, execution forbidden
	// tip  and  base  defined by Cell::Core, execution forbidden
	
	Mesh boundary ( );  // virtual from Cell::Core

	// two versions of  belongs_to  defined by Cell::Negative

	// glue_on_bdry_of  and  cut_from_bdry_of  defined by Cell::Core

	// methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'
	// see paragraph 10.1 in the manual

	// the two methods below are defined by Cell::Core, execution forbidden
	// void add_to_seg ( Cell::Positive::Segment * seg )
	// void remove_from_seg ( Cell::Positive::Segment * seg )

	void add_to_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_mesh ( Mesh::Core * msh ); // virtual from Cell::Core
  void add_to_bdry ( Mesh::Core * msh ); // virtual from Cell::Core
	void remove_from_bdry ( Mesh::Core * msh ); // virtual from Cell::Core

	// glue_on_my_bdry  and  cut_from_my_bdry  defined by Cell:Negative

#ifndef NDEBUG
	void print_everything ( ); // virtual from Cell::Core
#endif

}; // end of  class Cell::Negative::HighDim


//-----------------------------------------------------------------------------//
//------------------       core Meshes     ------------------------------------//
//-----------------------------------------------------------------------------//


class Mesh::Core

// represents a positive mesh
// negative meshes have no core (wrappers for negative meshes are built on-the-fly)

// abstract class, specialized in Mesh::ZeroDim, Mesh::Connected::{OneDim,HighDim},
// Mesh::MultiplyConnected::{OneDim,HighDim}, Mesh::Fuzzy, Mesh::STSI

{	public :

	Cell::Positive * cell_enclosed { nullptr };

	inline Core ( const tag::OfDimension &, size_t d , const tag::MinusOne & )  {	}
	
	virtual size_t get_dim_plus_one ( ) = 0;

	virtual size_t number_of ( const tag::Vertices & ) = 0;
	virtual size_t number_of ( const tag::Segments & ) = 0;
	virtual size_t number_of ( const tag::CellsOfDim &, size_t d ) = 0;

	// the four methods below are only relevant for connected one-dimensional meshes
	// so we forbid execution for now and then override them in Mesh::Connected::OneDim
	virtual Cell::Core * first_vertex ( );
	virtual Cell::Core * last_vertex ( );
	virtual Cell::Core * first_segment ( );
	virtual Cell::Core * last_segment ( );

	virtual void dispose () = 0;
	// remove_from this mesh each cell of maximum dimension
	// if the cell has no other meshes above, dispose of it
	// delete this (now empty) mesh
	
	//private :

	// the two methods below are only relevant for STSI meshes
	// so we forbid execution for now and then override them in Mesh::STSI
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                             ) const;

	virtual add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & ) = 0;
	virtual add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & ) = 0;
	virtual add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & ) = 0;
	virtual add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;
	virtual remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & ) = 0;

	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core *, size_t );
	// returns garbage; overriden by Mesh::Fuzzy and later by Mesh::STSI

	virtual void remove_from_cells
	( Cell::Core *, size_t, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// does nothing; overriden by Mesh::Fuzzy and later by Mesh::STSI
	
	// we are still in class Mesh::Core
	
	// iterators defined in iterator.cpp
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	inline CellIterator iterator
	( const tag::OverVertices, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
	{	return iterator ( tag::over_vertices, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_positive );  }
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & ) = 0;
	inline CellIterator iterator
	( const tag::OverVertices, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsNegative & )
	{	return iterator ( tag::over_vertices, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_negative );  }

	// we are still in class Mesh::Core
	
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	inline CellIterator iterator
	( const tag::OverSegments, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
	{	return iterator ( tag::over_segments, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_positive );  }
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & ) = 0;
	inline CellIterator iterator
	( const tag::OverSegments, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsNegative & )
	{	return iterator ( tag::over_segments, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_negative );  }
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsPositive & ) = 0;

	// we are still in class Mesh::Core
	
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive &             ) = 0;
	inline CellIterator iterator
	( const tag::OverSegments, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsPositive & )
	{	return iterator ( tag::over_cells_of_dim, d, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_positive );  }
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d,
	  const tag::ReverseOrder &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & ) = 0;
	virtual CellIterator iterator
	( const tag::OverCellsOfDim &, const size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative &             ) = 0;
	inline CellIterator iterator
	( const tag::OverSegments, const tag::ForcePositive &,
	  const tag::ReverseOrder &, const tag::ThisMeshIsNegative & )
	{	return iterator ( tag::over_cells_of_dim, d, tag::reverse_order,
                      tag::force_positive, tag::this_mesh_is_negative );  }
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsPositive & ) = 0;

#ifndef NDEBUG
	std::string name;
  virtual std::string get_name() = 0;
 	virtual void print_everything () = 0;
#endif

}; // end of  class Mesh::Core


class Mesh::ZeroDim : public Mesh::Core

// represents two points, one negative one positive

// zero-dimensional meshes only exist as boundary of a segment,
// which can be retrieved through the cell_enclosed attribute
// they are built on-the-fly

{	public :

	// attribute cell_enclosed inherited from Mesh::Core
	
	inline ZeroDim ( )
	:	Mesh::Core ( tag::of_dimension, 1, tag::minus_one )
	{ }

	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );
	// virtual from Mesh::Core, here execution forbidden
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	virtual void dispose ();  // virtual from Mesh::Core

	// private:

	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// the sixteen methods below are virtual from Mesh::Core, here execution forbidden
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	
	// we are still in class Mesh::ZeroDim
	
	// add_to_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	CellIterator iterator ( const tag::OverVertices &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first base (negative) then tip (positive)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first tip (positive) then base (negative)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first base then tip (both positive)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	// we iterate over the two vertices, first tip then base (both positive)

	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverVertices &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );

	CellIterator iterator ( const tag::OverVertices &, const tag::ThisMeshIsNegative & );
	// we iterate over the two vertices, first base (negative) then tip (positive)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );
	// we iterate over the two vertices, first tip (positive) then base (negative)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	// we iterate over the two vertices, first base then tip (both positive)

	CellIterator iterator
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	// we iterate over the two vertices, first tip then base (both positive)

	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverVertices &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );

	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverSegments &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverSegments &,, const tag::ForcePositive &,
	//    const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );

	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                   );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	  const tag::ThisMeshIsPositive &                                    );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive &       );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive &         );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ThisMeshIsNegative &                                   );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	  const tag::ThisMeshIsNegative &                                    );
	CellIterator iterator  // assert d == 0, call iterator over vertices
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative &       );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsNegative &         );

	// here is the core linking between cells and meshes
	// do not use directly
	// deep_connections_*** is called from add_to and remove_from
	static inline void make_deep_connections ( seg, Cell::Positive::Vertex * cll );
	static inline void break_deep_connections ( seg, Cell::Positive::Vertex * cll );
	static inline void make_deep_connections_rev ( seg, Cell::Positive::Vertex * cll );
	static inline void break_deep_connections_rev ( seg, Cell::Positive::Vertex * cll );

#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

};  // end of class Mesh::ZeroDim


class Mesh::Connected::OneDim : public Mesh::Core

// represents a connected positive mesh of dimension 1
// a chain of segments, either open or closed

{	public :

	// attribute cell_enclosed inherited from Mesh::Core

	size_t nb_of_segs;
	// useful for quickly answering to 'number_of'

	Cell::Negative::Vertex * first_ver;
	Cell::Positive::Vertex * last_ver;
	// if last_ver == first_ver->reverse, the mesh is a closed loop
	
	inline OneDim ( )
	:	Mesh::Core ( tag::of_dimension, 2, tag::minus_one )
	{ }

	OneDim ( const tag::Segment &, Cell::Negative::Vertex * A, Cell::Positive::Vertex * B,
	            const tag::DividedIn &, size_t n );
	// defined in global.cpp
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	Cell::Core * first_vertex ( );  // virtual from Mesh::Core, here overridden
	// returns a negative vertex
	Cell::Core * last_vertex ( );  // virtual from Mesh::Core, here overridden
	// returns a positive vertex
	Cell::Core * first_segment ( );  // virtual from Mesh::Core, here overridden
	Cell::Core * last_segment ( );  // virtual from Mesh::Core, here overridden

	virtual void dispose ();  // virtual from Mesh::Core

	// private:

	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage

	inline void order ( );
	// run over all segments, order them linearly, check that the mesh is connected
	// if it is a loop, set first_ver = Cell::ghost
	// if it is an open chain, set first_ver and last_ver
	// see paragraph 9.14 in the manual
	
#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

}; // end of  class Mesh::Connected::OneDim


class Mesh::Connected::HighDim : public Mesh::Core

// represents a connected positive mesh of dimension >= 2

{	public :

	// attribute cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	Cell start;
	
	inline HighDim ( const tag::OfDimension &, size_t dim_p1, const tag::MinusOne & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one )
	{	}

  HighDim ( const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA );
	// defined in global.cpp

	inline HighDim ( const tag::Quadrangle &,
		const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
		const tag::WithTriangles & wt = tag::not_with_triangles                       );

	inline HighDim ( const tag::Quadrangle &,
		const Cell & SW, const Cell & SE, const Cell & NE, const Cell & NW,
		size_t m, size_t n, const tag::WithTriangles & wt = tag::not_with_triangles );

	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	virtual void dispose ();  // virtual from Mesh::Core
	
	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage

	inline CellIterator iterator ( const tag::OverCellsOfDim &, size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::Connected::HighDim


class Mesh::MultiplyConnected::OneDim : public Mesh::Core

// represents a positive mesh of dimension 1 with several connected components

{	public :

	// attribute cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting/ending points

	inline OneDim ( )
	:	Mesh::Core ( tag::of_dimension, 2, tag::minus_one )
	{ }

  OneDim ( const tag::Segment &, Cell::Negative::Vertex * A,
		Cell::Positive::Vertex * B, const tag::DividedIn &, size_t n );
	// defined in global.cpp
	
	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	virtual void dispose ();  // virtual from Mesh::Core

	// private:
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage

	inline void order ( );
	// run over all segments, order them linearly, check that the mesh is connected
	// if it is a loop, set first_ver = Cell::ghost
	// if it is an open chain, set first_ver and last_ver
	// see paragraph 9.14 in the manual
	
#ifndef NDEBUG
	std::string get_name();  // virtual from Mesh::Core
	void print_everything ();  // virtual from Mesh::Core
#endif

}; // end of  class Mesh::MultiplyConnected::OneDim


class Mesh::MultiplyConnected::HighDim : public Mesh::Core

// represents a positive mesh of dimension >= 2 with several connected components

{	public :

	// attribute cell_enclosed inherited from Mesh::Core

	std::vector < size_t > nb_of_cells;
	// useful for quickly answering to 'number_of'

	// keep a list or a vector of starting cells
	
	inline HighDim ( const tag::OfDimension &, size_t dim_p1, const tag::MinusOne & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one )
	{	}

	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	virtual void dispose ();  // virtual from Mesh::Core
	
	static bool is_positive ( );
	static Mesh reverse ( Mesh::Core * core );

	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// add_to_cells ( Cell::Core *, size_t )  defined by Cell::Core, returns garbage

	inline CellIterator iterator ( const tag::OverCellsOfDim &, size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::MultiplyConnected::HighDim


class Mesh::Fuzzy : public Mesh::Core

// represents a positive mesh, disordered (includes 1D meshes)

// roughly speaking, a mesh is a collection of cells of the same dimension
// however, for efficiency purposes, we keep lists of cells of lower
// dimension as well; that's the purpose of the 'cells' vector :
// to keep lists of cells indexed by their dimension

{	public :
	
	// attribute cell_enclosed inherited from Mesh::Core

	// the 'cells' attribute holds lists of cells of 'this' mesh, indexed by their dimension
	// for maximum dimension, the cells are oriented
	// for lower dimension, the cells are always positive

	std::vector < std::list < Cell::Core* > > cells;

	inline Fuzzy ( const tag::OfDimension &, size_t dim_p1, const tag::MinusOne & )
	:	Mesh::Core ( tag::of_dimension, dim_p1, tag::minus_one ), cells ( d )
	{	}

	size_t get_dim_plus_one ( );  // virtual from Mesh::Core

	size_t number_of ( const tag::Vertices & );  // virtual from Mesh::Core
	size_t number_of ( const tag::Segments & );  // virtual from Mesh::Core
	size_t number_of ( const tag::CellsOfDim &, size_t d );  // virtual from Mesh::Core

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	void dispose ();  // virtual from Mesh::Core
	
	// private :
	
	// cell_in_front_of  and  cell_behind  defined by Mesh::Core, execution forbidden

	// the sixteen methods below are virtual from Mesh::Core
	// called from Cell::****tive::***::add_to_mesh and Cell::****tive::***::remove_from_mesh
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsNotBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsNotBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsNotBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsNotBdry & );
	void add_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void remove_pos_seg ( Cell::Positive::Segment *, const tag::MeshIsBdry & );
	void add_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void remove_neg_seg ( Cell::Negative::Segment *, const tag::MeshIsBdry & );
	void add_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void remove_pos_hd_cell ( Cell::Positive::HighDim *, const tag::MeshIsBdry & );
	void add_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );
	void remove_neg_hd_cell ( Cell::Negative::HighDim *, const tag::MeshIsBdry & );

	// add a cell to 'this->cells[d]' list, return iterator into that list
	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core * cll, size_t d ) override;
	
	// remove a cell from 'this->cells[d]' list using the provided iterator
	virtual void remove_from_cells
	( Cell::Core *, size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// virtual from Cell::Core, here overriden, later overriden again by Mesh::STSI
	
	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	CellIterator iterator
	( const tag::OverVertices &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverVertices &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator
	( const tag::OverVertices &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverVertices &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );
	CellIterator iterator
	( const tag::OverVertices &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverVertices &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverVertices &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator
	( const tag::OverSegments &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverSegments &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive & );
	CellIterator iterator
	( const tag::OverSegments &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );
	CellIterator iterator
	( const tag::OverSegments &, const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverSegments &, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative & );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverSegments &,, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsNegative & );
	CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsPositive & );
	CellIterator iterator  // execution forbidden
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ThisMeshIsPositive &                                   );
	CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	  const tag::ThisMeshIsPositive &                                    );
	CellIterator iterator  // execution forbidden
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsPositive &       );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsPositive &         );
	CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ThisMeshIsNegative & );
	CellIterator iterator  // execution forbidden
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ThisMeshIsNegative &                                   );
	CellIterator iterator
	( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	  const tag::ThisMeshIsNegative &                                    );
	CellIterator iterator  // execution forbidden
	( const tag::OverCellsOfDim &, size_t d, const tag::ReverseOrder &,
	  const tag::ForcePositive &, const tag::ThisMeshIsNegative &       );
	// CellIterator iterator  // inlined by Mesh::Core
	// ( const tag::OverCellsOfDim &, size_t d, const tag::ForcePositive &,
	//   const tag::ReverseOrder &, const tag::ThisMeshIsNegative &         );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::Fuzzy


class Mesh::STSI : public Mesh::Fuzzy

// represents a positive mesh, maybe self-touching, maybe self-intersecting
// (includes 1D meshes)

// mainly used for iterators over connected meshes
// also used for progressive mesh generation

{	public :
	
	// attribute cell_enclosed inherited from Mesh::Core

	// 'cells' attribute inherited from Mesh::Fuzzy

	// in 'singular' we keep pair of adjacent cells
	// the common face of such a pair is a singular face
	// that is, a face where the mesh touches itself
	std::vector < std::pair < Cell, Cell > > singular;
	
	inline STSI ( const tag::OfDimension &, size_t dim_p1, const tag::MinusOne & )
	:	Mesh::Fuzzy ( tag::of_dimension, dim_p1, tag::minus_one )
	{	}

	// size_t get_dim_plus_one ( )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Vertices & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::Segments & )  defined in Mesh::Fuzzy
	// size_t number_of ( const tag::CellsOfDim &, size_t d )  defined in Mesh::Fuzzy

	// first_vertex, last_vertex, first_segment and last_segment
	// are defined by Mesh::Core, execution forbidden

	void build_rectangle ( const Mesh & south, const Mesh & east,
		const Mesh & north, const Mesh & west, bool cut_rectangles_in_half );
	// defined in global.cpp

	virtual void dispose ();  // virtual from Mesh::Core
	
	static Mesh reverse ( Mesh::Core * core );

	// private :

	// the two methods below are defined in Cell::Core, we override them here
	virtual Cell::Core * cell_in_front_of
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;
	virtual Cell::Core * cell_behind
	( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
	  const tag::SurelyExists & se = tag::surely_exists                   ) const override;

	// add a cell to 'this->cells[d]' list, return iterator into that list
	virtual std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		add_to_cells ( Cell::Core * cll, size_t d ) override;
	// virtual from Cell::Core, overriden by Mesh::Fuzzy, here overriden again

	// remove a cell from 'this->cells[d]' list using the provided iterator
	virtual void remove_from_cells
	( Cell::Core *, size_t d, std::map<Mesh::Core*,Cell::field_to_meshes>::iterator );
	// virtual from Cell::Core, overriden by Mesh::Fuzzy, here overriden again

	// iterators are virtual from Mesh::Core and are defined in iterator.cpp

	inline CellIterator iterator ( const tag::OverCellsOfDim &, size_t d );
	
#ifndef NDEBUG
  std::string get_name();  // virtual from Mesh::Core
 	void print_everything ();  // virtual from Mesh::Core
#endif
	
}; // end of  class Mesh::STSI

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


// negative meshes aren't kept in the computer, they are just abstract concepts
// to be more precise: there are Mesh objects (wrappers) refering to negative meshes
// but their core points to a positive Mesh::Core; there are no negative Mesh::Cores
// the attribute Mesh::sign allows us to distinguish
// between a positive mesh and a negative one

class Mesh::Sign

{	public :
	virtual bool is_positive ( ) = 0;
	Class Positive;  Class Negative;                 };


class Mesh::Sign::Positive : public Mesh::Sign

{	public :
	bool is_positive ( );  // virtual from Mesh::Sign
};

class Mesh::Sign::Positive : public Mesh::Sign

{	public :
	bool is_positive ( );  // virtual from Mesh::Sign
};

//-----------------------------------------------------------------------------//
//-----------------------------------------------------------------------------//


inline Mesh::Mesh ( const tag::OfDimension &, size_t d, const tag::GreaterThanOne &,
                    const tag::IsPositive & ispos                                     )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	core { new Mesh::Positive ( tag::of_dimension, d+1, tag::minus_one ) },
	meth { & Mesh::Positive::methods_pos }
{	assert ( d > 1 );  }
	

inline Mesh::Mesh ( const tag::OfDimension &, size_t d, const tag::MightBeOne &,
                    const tag::IsPositive & ispos                                     )
// by default, ispos = tag::is_positive, so may be called with only three arguments
:	core { nullptr }, meth { & Mesh::Positive::methods_pos }
{	if ( d > 1 )
		this->core = new Mesh::Positive ( tag::of_dimension, d+1, tag::minus_one );
	else
	{	assert ( d == 1 );
		this->core = new Mesh::OneDim::Positive ( );  }                               }
	

inline Mesh::Mesh ( const tag::OfDimensionOne &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	core { new Mesh::OneDim::Positive ( ) }, meth { & Mesh::OneDim::Positive::methods_pos }
{ }
	

inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only two arguments
// used in Mesh::Negative::reverse and in Cell::boundary
: core { c }, meth { & Mesh::Positive::methods }
{	assert ( c );
	assert ( ( meth == & Mesh::Positive::methods_pos ) or
	         ( meth == & Mesh::OneDim::Positive::methods_pos ) );  }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::CellsSurelyExist &                                      )
// builds a negative mesh from a positive one, assuming that reverse cells exist
// used in Cell::boundary and in Mesh::Mesh below
: core { c }, meth { & Mesh::Negative::methods }
{	assert ( c );
#IFNDEF NDEBUG
	// check that all cells have reverse
#ENDIF
	assert ( ( meth == & Mesh::Positive::methods_neg ) or
	         ( meth == & Mesh::OneDim::Positive::methods_neg ) );  }


inline Mesh::Mesh ( const tag::WhoseCoreIs &, Mesh::Core * c, const tag::IsNegative &,
                    const tag::BuildCellsIfNec & b                                     )
// builds a negative mesh from a positive one, creating reverse cells if necessary
// used in Mesh::Positive::reverse
:	Mesh ( tag::whose_core_is, c, tag::is_negative, tag::cells_surely_exist )
{	std::list < Cell::Core * > & cls = c->cells.back();  // maximum dimension
	std::list<Cell::Core*>::iterator it = cls.begin(), it_e = cls.end();
	// cell[0] means cells of same dimension as the mesh
	for ( ; it != it_e; it++ )
	{	Cell::Core * cll_p = *it;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p );  assert ( cll_rev_p == cll_p->reverse_p );             }     }


inline Mesh::Mesh ( const tag::Segment &, const Cell & A, const Cell & B,
                    const tag::DividedIn &, size_t n                      )
: Mesh ( tag::whose_core_is, new Mesh::OneDim::Positive
  ( tag::segment, (Cell::Negative::Vertex*) A.core, (Cell::Positive::Vertex*) B.core,
    tag::divided_in, n ) )
{	}


inline Mesh::Mesh
( const tag::Pretty &, const tag::Segment &, const Cell & A, const Cell & B,
  const tag::DividedIn &, size_t n )
:	Mesh ( tag::of_dimension_one )
{	this->pretty_constructor ( tag::segment, A, B, tag::divided_in, n );  }


inline Mesh::Mesh
(	const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
: Mesh ( tag::whose_core_is, new Mesh::Positive ( tag::triangle, AB, BC, CA ) )
{	}


inline Mesh::Mesh
(	const tag::Pretty &, const tag::Triangle &, const Mesh & AB, const Mesh & BC, const Mesh & CA )
:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )
{	this->pretty_constructor ( tag::triangle, AB, BC, CA );  }


inline Mesh::Mesh
(	const tag::Pretty &, const tag::Quadrangle &,
	const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	const tag::WithTriangles & wt                                                  )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::of_dimension, 2, tag::greater_than_one )
{	this->pretty_constructor ( tag::quadrangle, south, east, north, west, wt );  }


inline Mesh::Mesh
(	const tag::Quadrangle &,
	const Mesh & south, const Mesh & east, const Mesh & north, const Mesh & west,
	const tag::WithTriangles & wt                                                  )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::whose_core_is, new Mesh::Positive
	        ( tag::quadrangle, south, east, north, west, wt ) )
{	}

inline Mesh::Mesh ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
                    const Cell & NE, const Cell & NW, size_t m, size_t n,
                    const tag::WithTriangles & wt                               )
// 'wt' defaults to 'tag::not_with_triangles',
:	Mesh ( tag::whose_core_is, new Mesh::Positive ( tag::quadrangle, SW, SE, NE, NW, m ,n, wt ) )
{	}


inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
  this->join_list ( l );                                         }

inline Mesh::Mesh
( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3, const Mesh & m4 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );
	l.push_back ( m3 );  l.push_back ( m4 );
  this->join_list ( l );                   }

inline Mesh::Mesh ( const tag::Join &, const Mesh & m1, const Mesh & m2, const Mesh & m3,
                                       const Mesh & m4, const Mesh & m5 )
:	Mesh ( tag::of_dim, m1.dim(), tag::might_be_one )
{	std::list < Mesh > l;
	l.push_back ( m1 );  l.push_back ( m2 );  l.push_back ( m3 );
	l.push_back ( m4 );  l.push_back ( m5 );
  this->join_list ( l );                                         }

inline Mesh::Mesh ( const tag::Join &, const std::list<Mesh> & l )
:	Mesh ( tag::of_dim, l.front().dim(), tag::might_be_one )
{	this->join_list ( l );  }


//-----------------------------------------------------------------------------//


inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh & msh )
:	Cell ( tag::whose_core_is, new Cell::Positive ( tag::whose_boundary_is, msh.core ) )
{	assert ( msh.is_positive() );  }


// inline Cell::Cell ( const tag::WhoseBoundaryIs &, Mesh::Core * msh )
// :	Cell ( tag::whose_core_is, new Cell::Positive ( tag::whose_boundary_is, msh ) )
// { }

inline Cell::Cell
( const tag::ReverseOf &, const Cell & direct_cell, const tag::BuildIfNotExists & build )
// 'build' defaults to tag::build_if_not_exists
// so constructor may be called with only two arguments
:	Cell ( tag::whose_core_is, direct_cell.core->reverse ( tag::build_if_not_exists ) )
{	}

inline Cell::Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::SurelyExists & )
:	Cell ( tag::whose_core_is, direct_cell.core->reverse_p )
{	}

inline Cell::Cell ( const tag::ReverseOf &, const Cell & direct_cell, const tag::MayNotExist & )
:	core ( direct_cell.core->reverse_p )
{	}


inline Cell::Cell ( const tag::BehindFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::SurelyExists & se )
// 'se' defaults to tag::surely_exists, so constructor may be called with only four arguments
:	core { msh.cell_behind ( face.core, tag::may_not_exist ) }
{	assert ( this->core );   }

inline Cell::Cell ( const tag::BehindFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::MayNotExist & )
:	core { msh.cell_behind ( face.core, tag::may_not_exist ) }
{	}

inline Cell::Cell ( const tag::InFrontOfFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::SurelyExists & se )
// 'se' defaults to tag::surely_exists, so constructor may be called with only four arguments
:	core { msh.cell_in_front_of ( face.core, tag::may_not_exist ) }
{	assert ( this->core );   }

inline Cell::Cell ( const tag::InFrontOfFace &, const Cell & face,
	const tag::WithinMesh &, const Mesh & msh, const tag::MayNotExist & )
:	core { msh.cell_in_front_of ( face.core, tag::may_not_exist ) }
{	}
		

inline Cell::Cell ( const tag::Vertex &, const tag::IsPositive & ispos )
// by default, ispos = tag::is_positive, so may be called with only one argument
:	Cell ( tag::whose_core_is, new Cell::Positive::Vertex )
{	}


inline Cell::Cell ( const tag::Segment &, const Cell & A, const Cell & B )
: Cell ( tag::whose_core_is, new Cell::Positive::Segment
( (Cell::Negative::Vertex*) A.core, (Cell::Positive::Vertex*) B.core ) )

{	assert ( not A.is_positive() );
	assert ( B.is_positive() );       }


inline Cell::Cell ( const tag::Triangle &, const Cell & AB, const Cell & BC, const Cell & CA )
:	Cell ( tag::whose_core_is, new Cell::Positive
            ( tag::triangle, AB.core, BC.core, CA.core ) )
{	}

inline Cell::Cell ( const tag::Quadrangle &, const Cell & AB, const Cell & BC,
                                             const Cell & CD, const Cell & DA )
:	Cell ( tag::whose_core_is, new Cell::Positive
            ( tag::quadrangle, AB.core, BC.core, CD.core, DA.core ) )
{	}

//-----------------------------------------------------------------------------//


inline size_t Mesh::number_of ( const tag::CellsOfDim &, size_t d ) const
{	return this->core->number_of ( tag::cells_of_dim, d );  }

inline size_t Mesh::number_of ( const tag::Segments & ) const
{	return this->number_of ( tag::cells_of_dim, 1 );  }

inline size_t Mesh::number_of ( const tag::Vertices & ) const
{	return this->number_of ( tag::cells_of_dim, 0 );  }


inline Cell Mesh::first_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_vertex() );
	else
	{	Cell::Core * tmp = this->core->last_vertex()->reverse_p;
		assert ( tmp );
		return Cell ( tag::whose_core_is, tmp );                  }      }

inline Cell Mesh::last_vertex ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_vertex() );
	else
	{	Cell::Core * tmp = this->core->first_vertex()->reverse_p;
		assert ( tmp );
		return Cell ( tag::whose_core_is, tmp );                   }     }

inline Cell Mesh::first_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->first_segment() );
	else
		return Cell ( tag::whose_core_is, this->core->last_segment()->reverse_p );  }

inline Cell Mesh::last_segment ( ) const
{	assert ( this->dim() == 1 );
	if ( this->is_positive() )
		return Cell ( tag::whose_core_is, this->core->last_segment() );
	else
		return Cell ( tag::whose_core_is, this->core->first_segment()->reverse_p );  }


inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::SurelyExists & se ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
{	return Cell ( tag::in_front_of_face, face, tag::within_mesh, *this, tag::surely_exists );  }

inline Cell Mesh::cell_in_front_of ( const Cell face, const tag::MayNotExist & ) const
{	return Cell ( tag::in_front_of_face, face, tag::within_mesh, *this, tag::may_not_exist );  }

inline Cell Mesh::cell_behind ( const Cell face, const tag::SurelyExists & se ) const
// 'se' defaults to tag::surely_exists, so method may be called with only one argument
{	return Cell ( tag::behind_face, face, tag::within_mesh, *this, tag::surely_exists );  }

inline Cell Mesh::cell_behind ( const Cell face, const tag::MayNotExist & ) const
{	return Cell ( tag::behind_face, face, tag::within_mesh, *this, tag::may_not_exist );  }


inline bool Cell::belongs_to ( const Mesh & msh, const tag::Oriented & ) const
// what about negative meshes ? haa ?!!	
{	return this->core->belongs_to ( msh.core, tag::oriented );  }

inline bool Cell::belongs_to ( const Mesh & msh, const tag::NotOriented & ) const
// what about negative meshes ? haa ?!!	
{	return this->core->belongs_to ( msh.core, tag::not_oriented );  }

//-----------------------------------------------------------------------------//


inline Mesh::Positive::Positive ( const tag::Quadrangle &, const Mesh & south,
  const Mesh & east, const Mesh & north, const Mesh & west, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

:	Mesh::Core ( tag::of_dimension, 3, tag::minus_one )
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	this->build_rectangle ( south, east, north, west, cut_rectangles_in_half );  }


inline Mesh::Positive::Positive ( const tag::Quadrangle &, const Cell & SW, const Cell & SE,
	const Cell & NE, const Cell & NW, size_t m, size_t n, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

:	Mesh::Core ( tag::of_dimension, 3, tag::minus_one )
	
{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	Mesh south ( tag::segment, SW.reverse(), SE, tag::divided_in, m );
	Mesh east  ( tag::segment, SE.reverse(), NE, tag::divided_in, n );
	Mesh north ( tag::segment, NE.reverse(), NW, tag::divided_in, m );
	Mesh west  ( tag::segment, NW.reverse(), SW, tag::divided_in, n );

	// when these four meshes go out of scope, their core should be disposed of

	this->build_rectangle ( south, east, north, west, cut_rectangles_in_half );  }

//-----------------------------------------------------------------------------//


inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::MayNotExist & ) const

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell::Core * face_rev = face_p->reverse_p;
	if ( face_rev == nullptr ) return nullptr;
	else return this->cell_behind ( face_rev, tag::may_not_exist );  }
	

inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell::Core * face_rev = face_p->reverse_p;
	assert ( face_rev );
	return this->cell_behind ( face_rev, tag::surely_exists );  }
	

inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::MayNotExist & ) const

// return the cell to which 'face' belongs, nullptr if we are facing the boundary

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_p->cell_behind_within.find ( this->core );
		if ( it == face_p->cell_behind_within.end() ) return nullptr;
			// nothing behind us, we are touching the boundary
		assert ( it->second );
		return it->second;                                                }
		// face_p->cell_behind_within[this->core]
	else
	{	Cell::Core * face_rev_p = face_p->reverse_p;
		assert ( face_rev_p );  // we are in a negative mesh, all faces must have reverse
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_rev_p->cell_behind_within.find ( this->core );
		if ( it == face_rev_p->cell_behind_within.end() )
			return nullptr;  // we are facing the boundary
		Cell::Core * cll_rev_p = it->second;
		assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
		return cll_rev_p->reverse_p;                                        }  }


inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SurelyExists & se ) const

// 'se' defaults to tag::surely_exists, so method may be called with only one argument

// return the cell to which 'face' belongs

{	assert ( this->dim() == face_p->get_dim() + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_p->cell_behind_within.find ( this->core );
		assert ( it != face_p->cell_behind_within.end() );
		assert ( it->second );
		return it->second;                                                }
		// face_p->cell_behind_within[this->core]
	else
	{	Cell::Core * face_rev_p = face_p->reverse_p;
		assert ( face_rev_p );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = face_rev_p->cell_behind_within.find ( this->core );
		assert ( it != face_rev_p->cell_behind_within.end() );
		Cell::Core * cll_rev_p = it->second;
		assert ( cll_rev_p );  assert ( cll_rev_p->reverse_p );
		return cll_rev_p->reverse_p;                                        }  }

// use :
// Cell f (...);
// Cell A ( tag::in_front_of_face, f, tag::within_mesh, msh );
// Cell::Core * A_p = msh.cell_in_front_of ( f.core );
// Cell::Core * A_p = msh.cell_in_front_of ( f.core, tag::may_not_exist );


// the four methods below are only relevant for STSI meshes

inline Cell::Core * Mesh::cell__in_front_of
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::MayNotExist &                                                       ) const
// return the cell towards which 'face' points, as seen from neighbour
{	return this->core->cell__in_front_of
		( face_p, tag::seen_from, neighbour, tag::may_not_exist );  }

inline Cell::Core * Mesh::cell_in_front_of
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se = tag::surely_exists                             ) const;
// return the cell towards which 'face' points, as seen from neighbour
{	return this->core->cell__in_front_of ( face_p, tag::seen_from, neighbour, se );  }

inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::MayNotExist &                                                       ) const
// return the cell to which 'face' belongs, as seen from neighbour
{	return this->core->cell_behind ( face_p, tag::seen_from, neighbour, tag::may_not_exist );  }

inline Cell::Core * Mesh::cell_behind
( const Cell::Core * face_p, const tag::SeenFrom &, const Cell::Core neighbour,
  const tag::SurelyExists & se = tag::surely_exists                             ) const;
// return the cell to which 'face' belongs, as seen from neighbour
{	return this->core->cell_behind ( face_p, tag::seen_from, neighbour, se );  }


#ifndef NDEBUG
inline void Mesh::print_everything ( )
{	if ( not is_positive() ) std::cout << "(negative Mesh) ";
	core->print_everything ();                                }
#endif

//-----------------------------------------------------------------------------//


inline bool Cell::is_positive ( ) const
{	return this->core->is_positive ( );  }

inline bool Mesh::is_positive ( ) const
{	return this->sign->is_positive ( );  }

inline Cell Cell::get_positive ( )
{	return Cell ( tag::whose_core_is, this->core->get_positive() );  }


inline size_t Cell::dim ( ) const
{	return this->core->get_dim ( );  }

inline size_t Mesh::dim ( ) const
{	return Mesh::diff ( this->core->get_dim_plus_one(), 1 );      }
// Mesh::diff  provides a safe way to substract two size_t numbers


inline Cell Cell::reverse ( const tag::BuildIfNotExists & build ) const
// 'build' defaults to tag::build_if_not_exists, so method may be called with no arguments
{	return Cell ( tag::reverse_of, *this, tag::build_if_not_exists );  }

inline Cell Cell::reverse ( const tag::MayNotExist & ) const
{	return Cell ( tag::reverse_of, *this, tag::may_not_exist );  }

inline Cell Cell::reverse ( const tag::SurelyExists & ) const
{	return Cell ( tag::reverse_of, *this, tag::surely_exists );  }

inline Mesh Mesh::reverse ( ) const
{	if ( this->is_positive() )
		return Mesh ( tag::whose_core_is, core, tag::is_negative, tag::build_cells_if_necessary );
	// else
	return Mesh ( tag::whose_core_is, core, tag::is_positive );                                   }
	

inline bool Cell::has_reverse ( ) const
{	return this->core->reverse_p;  }


inline Mesh Cell::boundary ( ) const

{	assert ( this->core );
	return this->core->boundary();    }

inline Cell Cell::tip () const
{	// assert ( this->core->tip() );
	return Cell ( tag::whose_core_is, this->core->tip() );  }

inline Cell Cell::base () const
{	// assert ( this->core->base() );
	return Cell ( tag::whose_core_is, this->core->base() );  }


#ifndef NDEBUG
inline void Cell::print_everything ( )
{	core->print_everything ( );  }
#endif


inline void Cell::glue_on_bdry_of ( Cell & cll )

// glue 'this' face on the boundary of cell 'cll'
// any of them may be negative

{	this->core->glue_on_bdry_of ( cll.core );  }


inline void Cell::cut_from_bdry_of ( Cell & cll )

// cut 'this' face from the boundary of cell 'cll'
// any of them may be negative

{	this->core->cut_from_bdry_of ( cll.core );  }


inline void Cell::add_to ( Mesh & msh )

// add 'this' cell to the mesh 'msh' by calling the virtual method add_to_mesh
// if 'msh' is the boundary of a cell, use instead 'glue_on_bdry_of'

{	assert ( this->dim() == msh.dim() );
	assert ( this->dim() > 0 );
	if ( msh.is_positive() )  this->core->add_to_mesh ( msh.core );
	else
	{	assert ( this->core->reverse_p );
		this->core->reverse_p->add_to_mesh ( msh.core );  }             }
// for negative Meshes, the core points towards the reverse Mesh::Core


inline void Cell::remove_from ( Mesh & msh )

// remove 'this' cell from the mesh 'msh' by calling the virtual method remove_from_mesh
// if 'msh' is the boundary of a cell, 'cut_from_bdry_of' should be used instead

{	assert ( this->dim() == msh.dim() );
	assert ( this->dim() > 0 );
	if ( msh.is_positive() )  this->core->remove_from_mesh ( msh.core );
	else
	{	assert ( this->core->reverse_p );
		this->core->reverse_p->remove_from_mesh ( msh.core );  }             }
// for negative Meshes, the core points towards the reverse Mesh::Core


//-----------------------------------------------------------------------------//


inline Cell::Positive::Vertex::Vertex ( )
: Cell::Positive ( tag::of_dim, 0, tag::size_meshes, Mesh::maximum_dimension_plus_one )
{	}


inline Cell::Negative::Vertex::Vertex
( const tag::ReverseOf &, Cell::Positive::Vertex * direct_ver_p )
: Cell::Negative ( tag::of_dim, 0, tag::reverse_of, direct_ver_p )
{	}


inline Cell::Positive::Segment::Segment
( Cell::Negative::Vertex * Aa, Cell::Positive::Vertex * Bb )

: Cell::Positive
		( tag::of_dim, 1, tag::size_meshes,
		  Util::assert_diff ( Mesh::maximum_dimension_plus_one, 1 ) ),
	// Util::assert_diff provides a safe way to substract two 'size_t' numbers
	base_p { Aa }, tip_p { Bb }

{	// below is a much simplified version of Negative::Vertex::add_to
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the base has already been correctly initialized
	assert ( Aa->reverse_p );
	Cell::Positive::Vertex * pos_Aa = Util::assert_cast
		< Cell::Core*, Cell::Positive::Vertex* > ( Aa->reverse_p );
	assert ( pos_Aa->meshes_same_dim.find(this) == pos_Aa->meshes_same_dim.end() );
	// pos_Aa->meshes[0][msh] = Cell::field_to_meshes { 0, 1 };
	// the third component 'where' is irrelevant here
	pos_Aa->meshes_same_dim.emplace ( std::piecewise_construct,
	      std::forward_as_tuple(msh), std::forward_as_tuple(-1) );
	// below is a much simplified version of Positive::Vertex::add_to
	// that's because 'this' segment has just been created, so it has no meshes above
	// also, the tip has already been correctly initialized
	assert ( Bb->meshes_same_dim.find(this) == Bb->meshes_same_dim.end() );
	// Bb->meshes[0][msh] = Cell::field_to_meshes { 1, 0 };
	// the third component 'where' is irrelevant here
	Bb->meshes_same_dim.emplace ( std::piecewise_construct,
	   std::forward_as_tuple(this), std::forward_as_tuple(1) );                         }


inline Cell::Negative::Segment::Segment
( const tag::ReverseOf &, Cell::Positive::Segment * direct_seg_p )

: Cell::Negative ( tag::of_dim, 1, tag::reverse_of, direct_seg_p )

// we must make sure that both extremities of 'direct_seg_p' have a reverse
// well, the base surely has one since it's a Negative::Vertex

{	assert ( direct_seg_p->base_p );
	assert ( direct_seg_p->tip_p );
	assert ( direct_seg_p->base_p->reverse_p );
	direct_seg_p->tip_p->reverse ( tag::build_if_not_exists );
	assert ( direct_seg_p->tip_p->reverse_p );                  }

		
inline Cell::Positive::HighDim::HighDim
( const tag::OfDimension &, size_t d, const tag::WhoseBoundaryIs &, Mesh::Core * msh )

:	Cell::Positive ( tag::of_dim, d, tag::size_meshes,
	                    Util::assert_diff ( Mesh::maximum_dimension_plus_one, d ) ),
	// Util::assert_diff provides a safe way to substract two 'size_t' numbers
	boundary_p ( msh )

{	assert ( msh );
	assert ( msh->get_dim_plus_one() == d );
	msh->cell_enclosed = this;                }


inline Cell::Positive::HighDim::HighDim
( const tag::WhoseBoundaryIs &, Mesh::Core * msh )
:	Cell::Positive ( tag::of_dim, msh->get_dim_plus_one(), tag::whose_bdry_is, msh )
{	}


inline Cell::Negative::HighDim::HighDim
( const tag::OfDimension, size_t d, const tag::ReverseOf &, Cell::Positive * direct_cell_p )
	
: Cell::Negative ( tag::of_dim, d, tag::reverse_of, direct_cell_p )

// we must make sure that all faces of 'direct_cell_p' have a reverse

{	assert ( direct_cell_p );
	assert ( direct_cell_p->get_dim() == d );
	assert ( direct_cell_p->boundary_p );
	std::list < Cell::Core * > & cls = direct_cell_p->boundary_p->cells.back();
	// cells of same dimension as the mesh (in this case, faces of 'this')
	std::list<Cell::Core*>::iterator it = cls.begin(), it_e = cls.end();
	for ( ; it != it_e; it++ )
	{	Cell::Core * cll_p = *it;
		assert ( cll_p );
		Cell::Core * cll_rev_p = cll_p->reverse ( tag::build_if_not_exists );
		assert ( cll_rev_p == cll_p->reverse_p );  assert ( cll_rev_p );             }    }


inline Cell::Negative::HighDim::HighDim
( const tag::ReverseOf &, Cell::Positive * direct_cell_p )
: Cell::Negative ( tag::of_dim, direct_cell_p->get_dim(), tag::reverse_of, direct_cell_p )
{	}


inline Cell::Positive::HighDim::HighDim ( const tag::Triangle &,
	Cell::Core * AB, Cell::Core * BC, Cell::Core * CA )
	
: Cell::Positive ( tag::whose_boundary_is, new Mesh::OneDim::Positive )

{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CA->base()->reverse_p );
	assert ( CA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to ( this->boundary_p );
	BC->add_to ( this->boundary_p );
	CA->add_to ( this->boundary_p );                }
	

inline Cell::Positive::HighDim::HighDim ( const tag::Quadrangle &,
	Cell::Core * AB, Cell::Core * BC, Cell::Core * CD, Cell::Core * DA )
	
: Cell::Positive ( tag::whose_boundary_is, new Mesh::OneDim::Positive )

{	assert ( AB->get_dim() == 1 );
	assert ( BC->get_dim() == 1 );
	assert ( CD->get_dim() == 1 );
	assert ( DA->get_dim() == 1 );
	assert ( AB->tip() == BC->base()->reverse_p );
	assert ( BC->tip() == CD->base()->reverse_p );
	assert ( CD->tip() == DA->base()->reverse_p );
	assert ( DA->tip() == AB->base()->reverse_p );
	// no need for glue_on_bdry_of : 'this' cell has just been created, it has no meshes above
	assert  ( this->boundary_p );
	AB->add_to ( this->boundary_p );
	BC->add_to ( this->boundary_p );
	CD->add_to ( this->boundary_p );
	DA->add_to ( this->boundary_p );                }


inline void Cell::Positive::glue_common ( Cell::Core * face )
	
{	if ( this->meshes.size() == 0 ) return;  // this never happens !!
	std::map < Mesh::Core*, Cell::field_to_meshes > & tm0 = this->meshes[0];
	// use meshes_same_dim !!
	// '0' means the same dimension as the cell
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = tm0.begin(), it_e = tm0.end();
	for ( ; it != it_e; ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
			face->cell_behind_within[msh] = this;
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			rev_face->cell_behind_within[msh] = other_cell;   }           }           }
	

inline void Cell::Core::cut_from_bdry_of ( Cell::Core * cll )
{	cll->cut_from_my_bdry ( this );   }


inline void Cell::Positive::cut_common ( Cell::Core * face )
	
{	if ( this->meshes.size() == 0 ) return;
	std::map < Mesh::Core *, Cell::field_to_meshes > & tm0 = this->meshes[0];
	// '0' means the same dimension as the cell
	std::map<Mesh::Core*,Cell::field_to_meshes>::iterator
		it = tm0.begin(), it_e = tm0.end();
	for ( ; it != it_e; ++it )
	{	Mesh::Core * msh = it->first;
		std::list<Cell::Core*>::iterator wh = it->second.where;
		Cell::Core * other_cell = *wh;  assert ( other_cell );
		if ( other_cell == this )  // orientations match
		{	assert ( face->cell_behind_within.find(msh) !=
		           face->cell_behind_within.end()        );
			assert ( face->cell_behind_within[msh] == this );
			face->cell_behind_within.erase(msh);                  }
		else  // mismatched orientations
		{	assert ( other_cell == this->reverse_p );
			Cell::Core * rev_face { face->reverse_p };
			assert ( rev_face );
			assert ( rev_face->cell_behind_within.find(msh) !=
		           rev_face->cell_behind_within.end()        );
			assert ( rev_face->cell_behind_within[msh] == other_cell );
			rev_face->cell_behind_within.erase(msh);                     }  }                  }
	

inline void Cell::Core::glue_on_bdry_of ( Cell::Core * cll )
{	cll->glue_on_my_bdry ( this );   }

//-----------------------------------------------------------------------------//


inline void Mesh::OneDim::Positive::order ( )

// run over all segments, order them linearly, check that the mesh is connected
// if it is a loop, set first_ver = Cell::ghost
// if it is an open chain, set first_ver and last_ver
// see paragraph 9.14 in the manual

{	if ( this->first_ver ) return;
	// if first_ver is not null, the mesh is already ordered

	std::list<Cell::Core*>::iterator it0 = this->cells[1].begin();  // will change !
	if ( it0 == this->cells[1].end() )  // empty mesh
	{	assert ( this->cells[1].size() == 0 );  // will change !
		this->first_ver = ( Cell::Negative::Vertex * ) Cell::ghost;
		return;                                                       }
	Cell::Core * seg = *it0;
	assert ( seg );
	assert ( seg->get_dim() == 1 );
	size_t counter = 1;
	Cell::Core * ver, * neg_ver;
	Cell::Core * seg_ini = seg;
	while ( true )
	{	ver = seg->tip();
		assert ( ver );
		assert ( ver->is_positive() );
		neg_ver = ver->reverse_p;
		if ( neg_ver == nullptr )  goto open_chain;
		assert ( not neg_ver->is_positive() );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = neg_ver->cell_behind_within.find ( this );
		if ( it == neg_ver->cell_behind_within.end() )  goto open_chain;
		seg = it->second;
		assert ( seg );
		if ( seg == seg_ini )  break;
		counter ++;                                              }
	// we assume there is no such loop having only one segment, OK ?
	if ( seg == seg_ini )  // we are dealing with a loop
	{	assert ( counter == this->cells[1].size() );  // may change !
		this->first_ver = ( Cell::Negative::Vertex * ) Cell::ghost;
		return;                                                      }
  open_chain :
	this->last_ver = ( Cell::Positive::Vertex* ) ver;
	seg = seg_ini;
	while ( true )
	{	neg_ver = seg->base();
		assert ( neg_ver );
		assert ( not neg_ver->is_positive() );
		ver = neg_ver->reverse_p;
		assert ( ver->is_positive() );
		std::map<Mesh::Core*,Cell::Core*>::const_iterator
			it = ver->cell_behind_within.find ( this );
		if ( it == ver->cell_behind_within.end() )  break;
		seg = it->second;
		assert ( seg );
		if ( seg == seg_ini )  break;
		counter ++;                                                       }
	assert ( counter == this->cells[1].size() );  // may change !
	this->first_ver = ( Cell::Negative::Vertex * ) neg_ver;                   }
	
//-----------------------------------------------------------------------------//


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



inline void Mesh::add_link_rev  // static
( Cell::Core * cll, Cell::Core * o_cll, Mesh::Core * msh, short int cp, short int cn )
// we just switch the two counters
{	Mesh::add_link ( cll, o_cll, msh, cn, cp );  }


inline void Mesh::remove_link  // static
( Cell::Positive * cll, Cell::Core * o_cll, Mesh::Core * msh, short int cp, short int cn )

{	assert ( cll );  assert ( o_cll );  assert ( msh );
	assert ( cll->is_positive() ); // assert ( msh->is_positive() );
	size_t cll_dim = cll->get_dim(),
	       msh_dim_p1 = msh->get_dim_plus_one();
	assert ( msh_dim_p1 > cll_dim );
	size_t dif_dim = msh_dim_p1 - cll_dim - 1;
	assert ( cll->meshes.size() > dif_dim );
	typedef std::map <Mesh::Core*, Cell::field_to_meshes> maptype;
	maptype & cmd = cll->meshes[dif_dim];
	maptype::iterator cmdfm = cmd.find(msh);
	assert ( cmdfm != cmd.end() );
	short int c_p = cmdfm->second.counter_pos -= cp;
	short int c_n = cmdfm->second.counter_neg -= cn;
	assert ( ( c_p >= 0 ) and ( c_n >= 0 ) );
	if ( ( c_p == 0 ) and ( c_n == 0 ) )
	{	std::list<Cell::Core*>::iterator w = cmdfm->second.where;
		msh->cells[cll_dim].erase(w);
		cmd.erase(cmdm);                                                 }
} // end of Mesh::remove_link


inline void Mesh::remove_link_rev  // static
( Cell::Core * cll, Cell::Core * o_cll,
  Mesh::Core * msh, short int cp, short int cn )
// we just switch the two counters
{	Mesh::remove_link ( cll, o_cll, msh, cn, cp );    }



}  // namespace maniFEM

#endif
// ifndef MANIFEM_MESH_H
