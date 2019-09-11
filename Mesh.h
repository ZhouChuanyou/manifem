// src/manifem/Mesh.h 2019.09.11

#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include "assert.h"
#include "Util.h"


namespace ManiFEM
{
	class Cell; class Segment; class Mesh; class AnyCell;
	class CellIterator; class MeshIterator;
	class Manifold; class ImplicitManifold;
	class NumericField; class BlockFieldBase; class OneDimField; class BlockField;
	class ComposedField; class MultiDimField; template <class VT> class PointerVector;
	class RelationFunc; class RelationFuncField; class RelationFuncConst;
	class FunctionOnAnyCell; class FieldOnAnyCell;
	class VariationalProblem; template <class U, class V>	class Equality;
	class FiniteElement; class Integrator;   }

namespace ManiFEM { namespace FunctionOnMesh {
		class baseFunction; class Function;
		class Constant; class Sum; class Product; class Fraction; class Power;
		class Unknown; class Test; class Mapping; class Composition;
		class Derivative; class Integral; class combinIntegrals; class Condition;
		extern Function zero, one, minus_one;
		Function & from_field ( NumericField &, std::string );
		Function & unknown ( NumericField &, std::string );
		Function & test ( Function & );
		inline Function & constant (double);
		inline baseFunction * sum ( baseFunction *, baseFunction * );
		inline baseFunction * sum ( baseFunction *, baseFunction *,
		                            baseFunction * );
		inline baseFunction * sum ( baseFunction *, baseFunction *,
		                            baseFunction *, baseFunction * );
		baseFunction * sum ( std::list < baseFunction * > );
		inline baseFunction * product ( baseFunction *, baseFunction * );
		inline baseFunction * product ( baseFunction *, baseFunction *,
		                            baseFunction * );
		inline baseFunction * product ( baseFunction *, baseFunction *,
		                            baseFunction *, baseFunction * );
		baseFunction * product ( std::list < baseFunction * > );
		inline baseFunction * fraction ( baseFunction *, baseFunction * );
		inline baseFunction * power ( baseFunction *, double );
		inline Mapping & mapping ( FunctionOnMesh::Function & );
		inline void prescribe_on ( Cell & cll );
}	}

namespace ManiFEM { namespace hidden {
		void for_hanging_nodes_2d ( Mesh & ambient_mesh, Cell & square,
			Cell & AB, Cell * & AB1_p, Cell * & AB2_p, Cell & center, double epsi );
	}	}

using namespace ManiFEM;


//////////////////////////////////////////////////////////////
//                                                          //
//          HHH    HHHHHHH   H       H        HHHH          //
//         H   H   H         H       H       H    H         //
//        H        H         H       H       H              //
//        H        HHHHH     H       H        H             //
//        H        H         H       H         HHH          //
//        H        H         H       H            H         //
//         H   H   H         H       H       H    H         //
//          HHH    HHHHHHH   HHHHHH  HHHHHH   HHHH          //
//                                                          //
//////////////////////////////////////////////////////////////


// A 'Cell' is a topological entity, like a finite element but without
// degrees of freedom or other stuff. The focus is on neighbourhood relations.

// A Cell of dimension zero is a point.
// A Cell of dimension one is a segment (see also the class Segment below).
// A Cell of dimension two may be a triangle, a quadrangle or some other polygon.
// A Cell of dimension three may be a tetrahedron, a cube or some other polyhedron.
// Cells of dimension four or higher may be constructed,
// but their usefulness is questionable.

// Cells may be positively or negatively oriented, see the attribute 'positive'.
// Which is which depends only on the construction process. The first one to be
// created is positive, the other one will be created when we call the 'reverse' method.

// A Cell is mainly defined by its boundary (which is a Mesh of lower dimension).
// The orientation of a Cell is nothing more than an orientation of its boundary;
// see the comments on orientation of Meshes below.


class ManiFEM::Cell

// most methods in this class are defined in Mesh.cpp

{	public:

	// attributes :
	// (perhaps they should be private, but I don't care much about that)

	//	static Cell * celula_especial; // !

	short int dim;   // (topological) dimension
	// normally it should be const, but since we want to re-use deleted cells
	// (keeping them in a pool) that's when the dimension may change (or not?)
	
	// hidden_boundary == NULL for points (cells of dimesion 0)
	Mesh* hidden_boundary;
	
	// When a (positive) cell is constructed, hidden_reverse == NULL.
	// When the 'reverse' method is invoked, a new (negative) cell is created,
	// and the attribute 'hidden_reverse' of both cells is set accordingly.
	Cell* hidden_reverse;

	// We give up the 'short int sign' attribute and replace it by a 'Cell* positive'.
	// If 'this == positive' then the cell is positive; otherwise, it is negative.
	Cell* positive;

	struct field_to_meshes
	{	short int counter_pos;
		short int counter_neg;
		std::list<Cell*>::iterator where;    };

	// The 'meshes' attribute keeps information about all meshes
	// "above" 'this' cell, that is, all meshes containing 'this'.
	// It is indexed over the dimension of the mesh and, for each mesh,
	// keeps a 'field_to_meshes' value, containing two counters
	// and an iterator into the 'cells' field of that particular mesh.
	// Of course this implies quite some amount of redundant information,
	// but this redundancy makes the classes fast, especially for remeshing.
	// To confirm : only positive meshes and positive cells are linked like that.
	//              a negative cell has empty 'meshes' map
	//              a negative mesh never appears in a 'meshes' map
	hidden::IndexedMap < Mesh*, field_to_meshes > meshes;

	// see the method 'Mesh::cell_in_front_of' for the significance of 'cell_behind_within'
	std::map < Mesh*, Cell* > cell_behind_within;

	// (we are still in class Cell)  we keep numeric values here :
	std::vector < double > real_heap;
	std::vector < size_t > size_t_heap;
	std::vector < short int > short_int_heap;

	static hidden::VectorOfIntegers < size_t >
		real_heap_size, size_t_heap_size, short_int_heap_size,
		real_heap_size_r, size_t_heap_size_r, short_int_heap_size_r ;
	// the first three are for positive cells, the other three for negative cells

	// 'hook' is a map where we keep miscelaneous (typeless) properties.
	// we make heavy use of explicit casting here
	std::map < std::string, void * > hook;
	// talvez mudar para um vector, e um mapa estático map <string, size_t> ?

	// a list of functions to be called each time a new cell is created
	// two lists, in fact, one for positive cells the other for negative cells
	static hidden::IndexedList < void(*)(Cell*,void*) > init_cell, init_cell_r;
	// more data can be passed to the above functions by using
	static hidden::IndexedList < void* > data_for_init, data_for_init_r;
	
#ifndef NDEBUG
	// the name only exists for debugging purposes
	std::string hidden_name;
#endif

	// constructors (do not use directly, use factory functions)
	
	Cell ( short int d, bool is_segment = false );
	// accepts d == 1 if second argument is true, but builds an inconsistent segment !
	// so do not use directly, use factory function Cell::segment instead
	// (see paragraph 5.6 of the manual)
	// bottom line : never use the second argument (only used once in Segment::Segment)

	Cell ( const tag::ReverseOf &, Cell & ); // for building the reverse cell
	// do not use directly, let 'reverse' do the job
	
	// factory functions :
	static Cell & point ();
	static Cell & segment ( Cell&, Cell& );
	static Cell & triangle ( Cell&, Cell&, Cell& );
	static Cell & rectangle ( Cell&, Cell&, Cell&, Cell& );
#ifndef NDEBUG
	// for cells with name :
	static Cell & point ( std::string );
	static Cell & segment ( std::string, Cell&, Cell& );
	static Cell & triangle ( std::string, Cell&, Cell&, Cell& );
	static Cell & rectangle ( std::string, Cell&, Cell&, Cell&, Cell& );
#endif
	
	// disallow copying :

	Cell ( const Cell & ) = delete;
	Cell& operator= ( const Cell & ) = delete;

	// we are still in class Cell

	// do not use for points, since they have no boundary
	inline Mesh & boundary ()
	{	assert ( hidden_boundary != NULL );
		return *hidden_boundary;               }
	
	inline Cell & reverse ( const tag::OnTheFly & otf = tag::on_the_fly );
	// usually called without any arguments
	// if 'hidden_reverse' is not NULL, return it (de-referenced)
	// if 'hidden_reverse' is NULL, a reverse cell is created "on the fly" and returned
	// defined after class ManiFEM::Segment

	// in the version below, 'reverse' returns a pointer which may be NULL
	inline Cell * reverse ( const tag::MayNotExist & ) const
	{	return hidden_reverse;  }
	
	inline bool is_positive () const {	return this == positive;  }

	inline bool belongs_to ( Mesh & msh ) const;

	// Method 'glue_on_bdry_of' is intensively used when building a mesh,
	// e.g. within factory functions in Cell class.
	// It glues 'this' cell to the boundary of 'cll'.
	void glue_on_bdry_of ( Cell & cll );
	// Method 'cut_from_bdry_of' does the reverse : cuts 'this' cell from
	// the boundary of 'cll'. Used mainly in remeshing.
	void cut_from_bdry_of ( Cell & cll );
	
	// Methods 'add_to' and 'remove_from' add/remove 'this' cell to/from the mesh 'msh'.
	// If 'msh' is the boundary of some cell, methods 'glue_on_bdry_of'
	// and 'cut_from_bdry_of' (see above) should be used instead.
	void add_to ( Mesh & msh );
	void remove_from ( Mesh & msh );

	// only for segments :
	inline Cell & tip ();
	inline Cell & base ();

#ifndef NDEBUG
	// 'name', 'print_meshes_above' and 'print_everything' only exist for debugging purposes
	std::string name ();
	void print_meshes_above ();
	void print_everything ();
#endif

	Cell & split ( const tag::InTwoSegments &, const tag::CellIsSegment & );
	// split 'this' segment in two smaller segments
	// after split, 'this' will be the first half of the segment
	// return the second half

	Cell & split ( const tag::InTwoTriangles &, const tag::At &,
	               Cell & ver, const tag::CellIsRectangle &       );
	// split 'this' rectangle along a diagonal
	// returns the newly created segment, beginning at 'ver'

	void split ( const tag::InFourRectangles &, const tag::CellIsRectangle &,
             const tag::WithIn &, Mesh & ambient_mesh, double epsi = 0.    );
	// split 'this' rectangle in four smaller rectangles, creating hanging nodes
	// degenerated triangles make the transition towards bigger rectangles

	void is_part_of ( Cell & cll );
	// state that 'this' cell must belong to all meshes where 'cll' belongs
	// used when we split a cell

	void discard ( bool first_time = true );
	// begin by discarding the boundary of 'this'
	// then remove 'this' from all meshes above it
	// finally, keep 'this' cell in a pool of cells for future use (to do)

	static Cell & cartesian_product_orient ( Cell & cell1, Cell & cell2,
		std::map < Cell*, std::map < Cell*, std::pair < Cell*, bool > > > & cartesian, bool revert );
	// do not use directly; use Mesh::cartesian_product

	inline MeshIterator iter_over ( const tag::Meshes &, const tag::Above &, const tag::OfSameDim &, const tag::Oriented & );
	// defined after class ManiFEM::MeshIterator
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &, const tag::OfDimPlusOne &, const tag::Oriented & );
	// defined after class ManiFEM::CellIterator
	
	// (we are still in class Cell)  data for iterators :
	// There are two kinds of iterators having 'Cell' as base, for now.
	// One of them is a MeshIterator; it iterates over all (oriented) meshes "above" a
	// given cell, having the same dimension as the cell.
	// The other one is a CellIterator; it iterates over all (oriented) cells "above" a
	// given cell, having dimension immediately higher.
	// They are almost identical, the second one returns the 'cell_enclosed' in the first one.
	struct data_for_iter_min_dim
	{	std::map <Mesh*, Cell::field_to_meshes> *mesh_map;
		std::map <Mesh*, Cell::field_to_meshes> :: iterator iter;  };
	static void reset_mesh_iter_min_dim ( MeshIterator* );
	static void advance_mesh_iter_min_dim ( MeshIterator* );
	static bool valid_mesh_iter_min_dim ( MeshIterator* );
	static Mesh * deref_iter_pos_min_dim ( MeshIterator* );
	static Mesh * deref_iter_neg_min_dim ( MeshIterator* );
	static void reset_cell_iter_min_dim ( CellIterator*, Cell* );
	static void advance_cell_iter_min_dim ( CellIterator *it );
	static bool valid_cell_iter_min_dim ( CellIterator *it );
	static Cell * deref_cell_iter_min_dim ( CellIterator *it );

}; // end of class Cell

//////////////////////////////////////////////////////////////////////////////////

class ManiFEM::Segment : public Cell

// A Segment is a Cell with two more attributes : hidden_base and hidden_tip

{	public :

	Cell* hidden_base;
	Cell* hidden_tip;

	// constructors
	inline Segment () : Cell ( 1, true ), hidden_base(NULL),
		hidden_tip(NULL) { }
	
	inline Segment ( const tag::ReverseOf &, const Segment & );
	// for building the reverse segment
	// do not use directly, let 'reverse' do the job

};

// we could have implemented 'base' and 'tip' as virtual methods
// but this would ocuppy more space in the memory, since each cell
// would have two pointers to these methods


inline Segment::Segment ( const tag::ReverseOf &, const Segment & direct_seg ) :

// for building the reverse segment

	Cell ( tag::reverse_of, * ( (Cell*) &direct_seg ) ),
	hidden_tip ( & ( direct_seg.hidden_base->reverse() ) ),
	hidden_base ( & ( direct_seg.hidden_tip->reverse() ) )

{	}


inline Cell & Cell::tip ()
{	assert ( dim == 1 );
	return * ( ((Segment*)this)->hidden_tip );  }

inline Cell & Cell::base ()
{	assert ( dim == 1 );
	return * ( ((Segment*)this)->hidden_base );  }


inline Cell & Cell::reverse ( const tag::OnTheFly & otf )

// usually called without any arguments ('otf' defaults to tag::on_the_fly)
// if the cell has already a reverse, return it
// otherwise, create on-the-fly the reverse cell, and return it

{	if ( this->hidden_reverse == NULL ) 
	{	// a negative cell has always a reverse,
		// so this must be a positive cell
		assert ( this->is_positive() );
		Cell *rev;
		if ( this->dim == 1 )
			rev = (Cell*) new Segment ( tag::reverse_of, * ( (Segment*) this ) );
		else rev = new Cell( tag::reverse_of, *this );
		return *rev;                                                             }
	// else :
	return * ( this->hidden_reverse );                                             }

//////////////////////////////////////////////////////////////////////////////////

class ManiFEM::AnyCell

// an abstract, generic cell which runs over a given mesh
// implicitly subject to a universal quantifier
// used for variational formulations
// see class RelationFuncField

{	public:
  
	size_t dim;

	// the 'environment' may be NULL, a pointer to a mesh or a pointer to a cell
	Mesh * env_mesh;
	Cell * env_cell;

	// constructors

	inline AnyCell ( Cell *cll )
	{	env_mesh = NULL;
		env_cell = cll;   }
	
	inline AnyCell ( Mesh *msh, size_t d )
	// msh may be NULL, which means "for all cells, in any mesh"
	{	env_mesh = msh;
		env_cell = NULL;
		dim = d;          }
  
};


/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//          H     H   HHHHHH    HHHH    H    H   HHHHHHH    HHHH           //
//          HH   HH   H        H    H   H    H   H         H    H          //
//          H H H H   H        H        H    H   H         H               //
//          H  H  H   HHHH      HHH     HHHHHH   HHHHH      HHH            //
//          H     H   H            H    H    H   H             H           //
//          H     H   H             H   H    H   H              H          //
//          H     H   H        H    H   H    H   H         H    H          //
//          H     H   HHHHHH    HHHH    H    H   HHHHHHH    HHHH           //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

// Roughly speaking, a mesh is a collection of cells of the same dimension.
// However, for efficiency purposes, we keep lists of cells of lower
// dimension as well. That's the purpose of 'IndexedList <Cell*> cells' :
// to keep lists of cells indexed by their dimension.
// This implies quite some amount of redundant information,
// but this redundancy makes the classes fast, especially for remeshing.

// An orientation of the mesh (described by its attribute 'sign')
// is nothing more than an orientation of each of its cells (of maximum dimension).
// But these orientations cannot be arbitrary, they must be compatible
// in the sense that a face common to two cells must be seen as positive
// from one of the cells and as negative from the other cell.
// As a consequence, for connected meshes there are only two possible orientations.
// Although nothing prevents a Mesh to be disconnected, we only allow for two orientations,
// a positive one (sign == 1) and a negative one (sign == -1).
// Which is which depends only on the construction process. The first one to be
// created is positive, the other one will be created when we call the 'reverse' method.
// In principle, negative meshes will appear only as boundaries of negative cells.


class ManiFEM::Mesh

{	public:

	// attributes :
	// (perhaps they should be private, but I don't care much about that)

	const short int dim;   // (topological) dimension

	// We give up the 'short int sign' attribute and replace it by 'Mesh * positive'.
	// If 'this == positive" then the cell is positive. Otherwise, it is negative
	Mesh* positive;
	
	// If the mesh is the boundary of a certain cell, [a pointer towards] that cell
	// is kept as 'cell_enclosed'. Otherwise, 'cell_enclosed' is NULL.
	Cell* cell_enclosed;
	
	// When a (positive) mesh is constructed, hidden_reverse == NULL.
	// When the 'reverse' method is invoked, a new (negative) mesh is created,
	// and the attribute 'hidden_reverse' of both meshes is set accordingly.
	Mesh* hidden_reverse;
	
	// The 'cells' attribute holds lists of cells of 'this' mesh,
	// indexed by their dimension. These are always positive cells.
	// The 'IndexedList' type is defined in the file 'IndexedPointer.h'.
	// Empty if 'this' mesh is negative.
	hidden::IndexedList < Cell* > cells;
	// podemos experimentar usar aqui uma                       //
	// map<Cell*,map<Mesh*,Cell:field_to_meshes>::iterator>     //
	// e ver se a velocidade aumenta significativamente         //

	// 'hook' is a map where we keep miscelaneous (typeless) properties.
	// we make heavy use of explicit casting here
	std::map <std::string, void*> hook;
	
	static short int intended_dimension;
	static Cell* current_cell_for_prescribe;
	static Manifold * environment;

	// constructors (do not use directly, use factory functions)
	
	inline Mesh ( short int d ) : dim (d), hidden_reverse (NULL),
		cell_enclosed (NULL), positive (this)
	{	 };

	Mesh ( const tag::ReverseOf &, Mesh & );  // for building the reverse mesh
	// defined in Mesh.cpp
	// do not use directly, let 'reverse' do the job

	// disallow copying :

	Mesh ( const Mesh & ) = delete;
	Mesh& operator= ( const Mesh & ) = delete;

	inline Mesh & reverse ( const tag::OnTheFly & otf = tag::on_the_fly );
	// usually called without any arguments
	// if 'hidden_reverse' is not NULL, return it
	// if 'hidden_reverse' is NULL, a reverse mesh is created "on the fly" and returned

	// in the version below, 'reverse' returns a pointer which may be NULL
	inline Mesh * reverse ( const tag::MayNotExist & ) const
	{	return hidden_reverse;  }
	
	inline bool is_positive () const {	return this == positive;  }

	inline size_t number_of ( const tag::Cells &, const tag::OfMaxDim & )
	{	return number_of ( tag::cells, tag::of_dim, dim );  }
	
	inline size_t number_of ( const tag::Cells &, const tag::OfDim &, short int d )
	{	if ( is_positive() ) return cells[d]->size();
		else return hidden_reverse->cells[d]->size();  }

	// (we are still in class Mesh)  factory functions :

	// 'build_boundary' builds (an empty mesh which will be) the boundary of the given cell
	inline static void build_boundary ( Cell& );

	// 'segment' builds a one-dimensional mesh with two extremities
	static Mesh & segment ( Cell&, Cell&, size_t );  // in global.cpp

	// 'rectangle' builds a two-dimensional mesh
	// arguments are : the four sides, in correct order and correctly oriented
	// the last argument chooses between rectangular cells and triangular ones
	// (triangles are obtained by cutting a rectangle in halves)
	static Mesh & rectangle  // defined in global.cpp
		( Mesh&, Mesh&, Mesh&, Mesh&, const tag::WithTriangles & wt = tag::not_with_triangles );
	
	// 'triangle' builds a two-dimensional mesh
	// arguments are : the three sides, in correct order and correctly oriented
	static Mesh & triangle ( Mesh&, Mesh&, Mesh& );  // defined in global.cpp
	
	// join several meshes :
	inline static Mesh & join ( Mesh & a, Mesh & b );
	inline static Mesh & join ( Mesh & a, Mesh & b, Mesh & c );
	inline static Mesh & join ( Mesh & a, Mesh & b, Mesh & c, Mesh & d );
	static Mesh & join ( const std::list<Mesh*> & l );
	// defined in global.cpp

	// 'f' is a face in 'this' mesh; f->dim == this->dim-1
	// Suppose 'f' belongs to the boundary of some cell 'cll'.
	// Then f->cell_behind_within[this] is equal to 'cll'.
	// If f has reverse, then f->hidden_reverse->cell_behind_within[this]
	// gives the other cell, neighbour to 'cll', having 'f' as common face.
	// This is exactly what cell_in_front_of implements.
	// Interpretation : the faces are looking outwards, to the neighbour cells.
	inline Cell * cell_in_front_of ( Cell &, const tag::MayNotExist & );
	inline Cell * cell_behind ( Cell &, const tag::MayNotExist & );
	// return nullpointer if there is no cell (we are on the boundary of the mesh)
	inline Cell & cell_in_front_of ( Cell &, const tag::SurelyExists & se = tag::surely_exists );
	inline Cell & cell_behind ( Cell&, const tag::SurelyExists & se = tag::surely_exists );

	// for a cell of maximum dimension,
	// returns a cell with orientation compatible with 'this' orientation
	inline Cell & orient_cell ( Cell & cll );

	// (we are still in class Mesh)  output :
	
	void export_msh ( std::string f, std::map<Cell*,size_t> & ver_numbering );
	void export_msh ( std::string f );
	void export_msh ( std::string f, FiniteElement & fe );
	void draw_ps ( std::string f );
#ifndef NDEBUG
	// 'name', 'print_cells' and 'print_everything' only exist for debugging purposes
	std::string name ();
	void print_cells ();
	void print_everything ();
#endif
	
	// methods action_* should not be used directly in programs
	// they are passed to 'deep_connections' from  Cell::add_to and Cell::remove_from
	static void action_add ( Cell&, Mesh&, short int, short int );
	static void action_remove ( Cell&, Mesh&, short int, short int );
	static void action_add_rev ( Cell & cll, Mesh & msh, short int, short int );
	static void action_remove_rev ( Cell & cll, Mesh & msh, short int, short int );

	// here is the core linking between cells and meshes
	// do not use directly; this is called from add_to and remove_from
	void deep_connections ( Cell&, void (*action) ( Cell&, Mesh&, short int, short int ) );

	static Mesh & cartesian_product ( Mesh & mesh1, Mesh & mesh2 );

	void discard ();
	// first discard all cells
	// then keep 'this' mesh in a pool of cells for future use (to do)

	// iterator over Cells belonging to 'this'
	inline CellIterator iter_over ( const tag::Cells &, const tag::OfMaxDim &, const tag::Oriented & );
	// code example :  auto it = msh.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
	//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }
	inline CellIterator iter_over
		( const tag::Cells &, const tag::OfDim &, const short int, const tag::NotOriented & );
	// code example :  auto it = msh.iter_over ( tag::cells, tag::of_dim, d, tag::not_oriented );
	//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusOne &, const tag::Along &             );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusOne &, const tag::Along &, const tag::Reverse & );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusOne &, const tag::Around &            );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusOne &, const tag::Around &, const tag::Reverse & );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusTwo &, const tag::Along &             );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusTwo &, const tag::Along &, const tag::Reverse & );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusTwo &, const tag::Around &            );
	inline CellIterator iter_over ( const tag::Cells &, const tag::Above &,
		Cell & cll, const tag::OfDimPlusTwo &, const tag::Around &, const tag::Reverse & );
	// code example :  auto it = msh.iter_over ( tag::cells, tag::above, cll,
	//                                           tag::of_dim_plus_two, tag::along );
	//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }
	// only for one-dimensional meshes :
	inline CellIterator iter_over ( const tag::Vertices &, const tag::Around & );
	inline CellIterator iter_over ( const tag::Vertices &, const tag::Along & );
	inline CellIterator iter_over ( const tag::Segments &, const tag::Around & );
	inline CellIterator iter_over ( const tag::Segments &, const tag::Along & );
	inline CellIterator iter_over ( const tag::Vertices &, const tag::Around &,
	                                const tag::Reverse &                         );
	inline CellIterator iter_over ( const tag::Vertices &, const tag::Along &,
	                                const tag::Reverse &                         );
	inline CellIterator iter_over ( const tag::Segments &, const tag::Around &,
	                                const tag::Reverse &                         );
	inline CellIterator iter_over ( const tag::Segments &, const tag::Along &,
	                                const tag::Reverse &                         );
	// code example :  auto it = msh.iter_over ( tag::segments_along, tag::reverse );
	//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }
	// defined after class ManiFEM::CellIterator
	
	// (we are still in class Mesh)  data for iterators :
	// there are several types of iterator having a 'Mesh' as base
	// a CellIterator which iterates over all (oriented) cells of a given mesh,
	// of maximum dimension
	struct data_for_iter_max_dim
	{	std::list <Cell*> *cell_list;
		std::list <Cell*> :: iterator iter;  };
	// another one which iterates over non-oriented cells of a given mesh, of given dimension :
	struct data_for_iter_given_dim
	{	short int dim;
		std::list <Cell*> *cell_list;
		std::list <Cell*> :: iterator iter;  };
	// contour iterators include "segments around", "vertices around",
	// "segments along" and "vertices along"
	struct data_for_iter_contour
	{	Cell * current_vertex; Mesh *oriented_base;
		Cell * current_seg, * first_seg, * center;
		bool first_try, seg_got_out;                 };
	static void reset_iter_max_dim ( CellIterator*, Cell* );
	static void advance_iter_max_dim ( CellIterator* );
	static bool valid_iter_max_dim ( CellIterator* );
	static Cell * deref_iter_pos_max_dim ( CellIterator* );
	static Cell * deref_iter_neg_max_dim ( CellIterator* );
	static void reset_iter_around ( CellIterator*, Cell* );
	static void reset_iter_around_rev ( CellIterator*, Cell* );
	static void reset_iter_along_ver ( CellIterator*, Cell* );
	static void reset_iter_along_ver_rev ( CellIterator*, Cell* );
	static void reset_iter_along_seg ( CellIterator*, Cell* );
	static void reset_iter_along_seg_rev ( CellIterator*, Cell* );
	static void advance_iter_contour ( CellIterator* );
	static void advance_iter_contour_rev ( CellIterator* );
	static Cell * deref_iter_segs_contour ( CellIterator* );
	static Cell * deref_iter_vertices_contour (CellIterator*);
	static bool valid_iter_around ( CellIterator* );
	static bool valid_iter_segs_along ( CellIterator* );
	static bool valid_iter_vertices_along ( CellIterator* );
	static void reset_iter_given_dim ( CellIterator*, Cell* );
	static void advance_iter_given_dim ( CellIterator* );
	static bool valid_iter_given_dim ( CellIterator* );
	static Cell * deref_iter_given_dim ( CellIterator* );
	static void reset_iter_around_dual ( CellIterator*, Cell* );
	static void reset_iter_along_dual ( CellIterator*, Cell* );
	static void reset_iter_along_dual_rev ( CellIterator*, Cell* );
	static void advance_iter_dual ( CellIterator* );
	static void advance_iter_dual_rev ( CellIterator* );
	static bool valid_iter_along_dual ( CellIterator* );
	static bool valid_iter_around_dual ( CellIterator* );
	inline Cell * get_tri_for_iter_dual ( Cell * center );

	static inline void prescribe_on ( Cell* cll, void* )
	{ FunctionOnMesh::prescribe_on ( *cll );  }

}; // end of class ManiFEM::Mesh


inline void Mesh::build_boundary ( Cell & cll )

// builds (an empty mesh which will be) the boundary of cll

{	Mesh *msh = new Mesh ( cll.dim - 1 );
	msh->cell_enclosed = & cll;
	cll.hidden_boundary = msh;             }


inline Mesh & Mesh::reverse ( const tag::OnTheFly & otf )

// usually called without any arguments ('otf' defaults to tag::on_the_fly)
// if 'this' mesh has already a reverse, return it
// otherwise, create on-the-fly the reverse mesh, and return it

{	if ( this->hidden_reverse == NULL )
	{	// a negative mesh has always a reverse,
		// so this must be a positive mesh
		assert ( this->is_positive() );
		return * ( new Mesh ( tag::reverse_of, *this ) );  }
	else return * ( this->hidden_reverse );                 }


inline Cell * Mesh::cell_behind ( Cell & face, const tag::MayNotExist & )

// return the cell behind 'face'
// recall that the faces of a cell are looking outwards

{	short int d = this->dim;
	assert ( d == face.dim + 1 );
	if ( this->is_positive() )
	{	std::map<Mesh*,Cell*>::iterator it = face.cell_behind_within.find(this);
		if ( it == face.cell_behind_within.end() ) return NULL;
		return it->second;                                                        }
		// face.cell_behind_within[this]
	else
	{	assert ( this->hidden_reverse != NULL );
		Cell *cll = this->hidden_reverse->cell_in_front_of ( face, tag::may_not_exist );
		if ( cll ) return cll->hidden_reverse;
		else return NULL;                                                                 }
} // end of *Mesh::cell_behind


inline Cell & Mesh::cell_behind ( Cell & face, const tag::SurelyExists & se )

// usually called with one argument only ('se' defaults to tag::surely_exists)

{	Cell * cll = this->cell_behind ( face, tag::may_not_exist );
	assert ( cll != NULL );
	return *cll;                                                    }


inline Cell * Mesh::cell_in_front_of ( Cell & face, const tag::MayNotExist & )

// return the cell towards which 'face' is looking
// recall that the faces of a cell are looking outwards

{	Cell * fr = face.hidden_reverse;
	if ( fr == NULL ) return NULL;
	else return this->cell_behind ( * fr, tag::may_not_exist );
} // end of Mesh::cell_in_front_of


inline Cell & Mesh::cell_in_front_of ( Cell & face, const tag::SurelyExists & se )

// usually called with one argument only ('se' defaults to tag::surely_exists)

{	Cell * cll = this->cell_in_front_of ( face, tag::may_not_exist );
	assert ( cll != NULL );
	return *cll;                                                         }


inline Cell & Mesh::orient_cell ( Cell & c )

// for a cell of maximum dimension inside 'this' mesh,
// returns a cell with orientation compatible with 'this' orientation

{	assert ( c.dim == this->dim );
	bool ok = true;
	Mesh * msh = this;
	if ( ! ( msh->is_positive() ) )
	{	ok = false; msh = msh->hidden_reverse;  }
	assert ( msh != NULL );
	assert ( msh->is_positive() );
	Cell *cll = &c;
	if ( ! ( cll->is_positive() ) ) cll = cll->hidden_reverse; // ok = not ok ?
	assert ( cll != NULL );
	assert ( cll->is_positive() );
	Cell::field_to_meshes & field = (*(cll->meshes[msh->dim]))[msh];
	if (field.counter_pos != 0)
	{	assert (field.counter_pos == 1);
		assert (field.counter_neg == 0);   }
	else
	{	assert (field.counter_pos == 0);
		assert (field.counter_neg == 1);
		assert (cll->hidden_reverse != NULL);
		ok = not ok;                            }
	if (ok) return * cll;
	assert (cll->hidden_reverse != NULL);
	return * ( cll->hidden_reverse );                                      }


inline Mesh & Mesh::join ( Mesh & a, Mesh & b ) // static
{	std::list < Mesh * > l {&a,&b};
	return Mesh::join (l);           }

inline Mesh & Mesh::join ( Mesh & a, Mesh & b, Mesh & c )  // static
{	std::list < Mesh * > l {&a,&b,&c};
	return Mesh::join (l);              }

inline Mesh & Mesh::join ( Mesh & a, Mesh & b, Mesh & c, Mesh & d )  // static
{	std::list < Mesh * > l {&a,&b,&c,&d};
	return Mesh::join (l);                 }


inline bool Cell::belongs_to ( Mesh & msh ) const

// important : if this->dim < msh.dim then the orientations do not matter
// (in this case, usually the cell appears more than once in the mesh,
//   with different orientations)
// but if this->dim == msh.dim then the orientation matters
	
{	short int d = msh.dim;
	Cell * positive_this = this->positive;
	Mesh * positive_msh = msh.positive;
	bool same_orien = ( ( this == positive_this ) == ( (&msh) == positive_msh ) );
	std::map < Mesh*, Cell::field_to_meshes > * ptmd = positive_this->meshes[d];
	std::map < Mesh*, Cell::field_to_meshes > :: iterator
		it = ptmd->find(positive_msh);
	if ( it == ptmd->end() ) return false;
	if ( this->dim < d ) // orientation does not matter
		return true;
	assert ( ( ( it->second.counter_pos == 1 ) and ( it->second.counter_neg == 0 ) ) or
	         ( ( it->second.counter_pos == 0 ) and ( it->second.counter_neg == 1 ) ) );
	return same_orien == ( it->second.counter_pos == 1 );                                 }
	

inline Cell * Mesh::get_tri_for_iter_dual ( Cell * center )

// returns a triangle having 'center' as vertex, within 'this' mesh
// well, not necessarily a triangle, just any two-dimensional cell

{	center = center->positive;
	assert ( this->dim == center->dim + 2 );
	std::map<Mesh*,Cell::field_to_meshes> * cmd = center->meshes[this->dim-1];
	std::map<Mesh*,Cell::field_to_meshes>::iterator
	  it = cmd->begin(),  it_e = cmd->end();
	for ( ; it != it_e; it++ )
	{	Mesh * tri_bdry = it->first;
		Cell * tri = tri_bdry->cell_enclosed;
		if ( tri == NULL ) continue;
		// we must check that tri actually belongs to 'this' mesh
		std::map<Mesh*,Cell::field_to_meshes> * tmd = tri->meshes[this->dim];
		if ( tmd->find(this) == tmd->end() ) continue;
		// 'tri' must be positive because 'tri_bdry' is
		return tri;                                                            }
	assert ( false );                                                          }


		
////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//      H  HHHHHHHHH  HHHHHHH   HHHH       H    HHHHHHHHH    HHH     HHHH      HHH        //
//      H      H      H         H   H      H        H       H   H    H   H    H   H       //
//      H      H      H         H   H     H H       H      H     H   H   H     H          //
//      H      H      HHHH      HHHH      H H       H      H     H   HHHH       H         //
//      H      H      H         H H      H   H      H      H     H   H H         H        //
//      H      H      H         H  H     HHHHH      H      H     H   H  H         H       //
//      H      H      H         H   H   H     H     H       H   H    H   H     H  H       //
//      H      H      HHHHHHH   H    H  H     H     H        HHH     H    H     HH        //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////


// Iterators are classified into two types.
// MeshIterators iterate over lists of meshes, e.g. over all meshes above a given cell.
// CellIterators iterate over lists of cells,
// e.g. over all cells of given dimension belonging to a given mesh,
// or over all cells of given dimension "above" another, given, cell.


class ManiFEM::MeshIterator

{	public :
	void *base;
	void *data;
	Mesh *value;

	// pointers to functions :
	void (*re_set) (MeshIterator*);
	void (*ad_vance) (MeshIterator*);
	bool (*va_lid) (MeshIterator*);
	Mesh * (*dereference) (MeshIterator*);

	// constructors
	
	MeshIterator () = delete;
	
	// iterates over Meshes "above" 'cll'
	inline MeshIterator
		( const Cell & cll, const tag::Meshes &, const tag::Above &, const tag::OfSameDim &, const tag::Oriented & )
		: base ( (void*) cll.positive )
	{ };
	
	inline void reset()          {  re_set   (this);  }
	inline void advance()        {  assert ( va_lid (this) );  ad_vance (this);  }
	inline void operator++()     {  assert ( va_lid (this) );  ad_vance (this);  }
	inline void operator++(int)  {  assert ( va_lid (this) );  ad_vance (this);  }
	inline Mesh & operator*()
	{	assert ( va_lid (this) );
		Mesh * msh = dereference (this);
		assert ( msh ); // msh not NULL
		return * msh;                     }
	inline bool in_range()       {  return va_lid (this);  }
	inline void skip_null();  // useless

}; // end of class ManiFEM::MeshIterator


inline MeshIterator Cell::iter_over ( const tag::Meshes &, const tag::Above &, const tag::OfSameDim &, const tag::Oriented & )

{	MeshIterator it ( *this, tag::meshes, tag::above, tag::of_same_dim, tag::oriented );
	it.data = (void*) new Cell::data_for_iter_min_dim;
	it.re_set = Cell::reset_mesh_iter_min_dim;
	it.ad_vance = Cell::advance_mesh_iter_min_dim;
	it.va_lid = Cell::valid_mesh_iter_min_dim;
	if ( this->is_positive() ) it.dereference = Cell::deref_iter_pos_min_dim;
	else
	{	assert ( this->hidden_reverse != NULL );
		it.dereference = Cell::deref_iter_neg_min_dim;    }
	return it;                                                                        }
	

inline void MeshIterator::skip_null ()  // should be eliminated

{	while ( this->in_range() )
	{	this->value = this->dereference (this);
		if ( this->value != NULL ) break;
		assert ( false );
		this->advance ();                         }  }

////////////////////////////////////////////////////////////////////////////////////////

class ManiFEM::CellIterator

{	public :
	
	void *base;
	void *data;

	// pointers to functions :
	void (*re_set) ( CellIterator*, Cell* );
	void (*ad_vance) ( CellIterator* );
	bool (*va_lid) ( CellIterator* );
	Cell * (*dereference) ( CellIterator* );
	
	// constructors (do not use directly)
	// use instead pseudo-constructor Mesh::iter_over
	
	CellIterator () = delete;
	
	// iterate over Cells belonging to 'msh'
	// used in  Mesh::iter_over ( 18 possibilities )
	inline CellIterator ( const Mesh & msh ) : base ( (void*) msh.positive )
	{ };

	// iterates over Cells "above" 'cll'
	inline CellIterator ( const Cell & cll ) : base ( (void*) cll.positive )
	// anyway, 'base' is irrelevant here (never used),
	// see Cell::iter_over ( tag::Cells &, const tag::Above &, const tag::OfDimPlusOne &, tag::Oriented & )
	{ };
	
	inline void reset ()  {  re_set ( this, NULL );  }
	inline void reset ( Cell & p )  {  re_set ( this, & p );  }
	inline void advance()        {  assert ( va_lid (this) );  ad_vance (this);  }
	inline void operator++()     {  assert ( va_lid (this) );  ad_vance (this);  }
	inline void operator++(int)  {  assert ( va_lid (this) );  ad_vance (this);  }
	inline Cell & operator*()
	{	assert ( va_lid (this) );
		Cell * cll = dereference (this);
		assert ( cll ); // cll not NULL
		return * cll;                     }
	inline bool in_range()       {  return va_lid (this);  }
	inline void skip_null();
	
}; // end of class ManiFEM::CellIterator



inline CellIterator Cell::iter_over ( const tag::Cells &, const tag::Above &, const tag::OfDimPlusOne &, const tag::Oriented & )

// iterator over Cells "above" 'this' cell

{	CellIterator it ( *this );
//	      ( this, "cells above, of minimum dimension, oriented") );  }
	// the only relevant datum is another iterator, over meshes "above" cll
	MeshIterator * it_m = new MeshIterator  // !
		( *this, tag::meshes, tag::above, tag::of_same_dim, tag::oriented );
	it_m->data = (void*) new Cell::data_for_iter_min_dim;
	it_m->re_set = Cell::reset_mesh_iter_min_dim;
	it_m->ad_vance = Cell::advance_mesh_iter_min_dim;
	it_m->va_lid = Cell::valid_mesh_iter_min_dim;
	if ( this->is_positive() ) it_m->dereference = Cell::deref_iter_pos_min_dim;
	else
	{	assert ( this->hidden_reverse != NULL );
		it_m->dereference = Cell::deref_iter_neg_min_dim;    }
	it.data = (void*) it_m;
	it.re_set = Cell::reset_cell_iter_min_dim;
	it.ad_vance = Cell::advance_cell_iter_min_dim;
	it.va_lid = Cell::valid_cell_iter_min_dim;
	it.dereference = Cell::deref_cell_iter_min_dim;
	return it;                                                                      }


inline CellIterator Mesh::iter_over	( const tag::Cells &, const tag::Above &,
	Cell & cll, const tag::OfDimPlusOne &, const tag::Along &    )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_dual;
	it.ad_vance = Mesh::advance_iter_dual;
	it.va_lid = Mesh::valid_iter_along_dual;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over ( const tag::Cells &, const tag::Above &, Cell & cll,
	  const tag::OfDimPlusOne &, const tag::Along &, const tag::Reverse &    )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_dual_rev;
	it.ad_vance = Mesh::advance_iter_dual_rev;
	it.va_lid = Mesh::valid_iter_along_dual;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll, const tag::OfDimPlusOne &, const tag::Around & )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_dual;
	it.ad_vance = Mesh::advance_iter_dual;
	it.va_lid = Mesh::valid_iter_around_dual;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll,
	  const tag::OfDimPlusOne &, const tag::Around &, const tag::Reverse & )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_dual;
	it.ad_vance = Mesh::advance_iter_dual_rev;
	it.va_lid = Mesh::valid_iter_around_dual;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll, const tag::OfDimPlusTwo &, const tag::Along & )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_dual;
	it.ad_vance = Mesh::advance_iter_dual;
	it.va_lid = Mesh::valid_iter_along_dual;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll,
	  const tag::OfDimPlusTwo &, const tag::Along &, const tag::Reverse & )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_dual_rev;
	it.ad_vance = Mesh::advance_iter_dual_rev;
	it.va_lid = Mesh::valid_iter_along_dual;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll, const tag::OfDimPlusTwo &, const tag::Around & )
		
{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_dual;
	it.ad_vance = Mesh::advance_iter_dual;
	it.va_lid = Mesh::valid_iter_around_dual;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::Above &, Cell & cll,
	  const tag::OfDimPlusTwo &, const tag::Around &, const tag::Reverse & )

{	CellIterator it ( *this );
	assert ( this->dim == cll.dim + 2 );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	d->center = &cll;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_dual;
	it.ad_vance = Mesh::advance_iter_dual_rev;
	it.va_lid = Mesh::valid_iter_around_dual;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over ( const tag::Cells &, const tag::OfMaxDim &, const tag::Oriented & )

// iterator over Cells belonging to 'this' mesh
// code example :  auto it = msh.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	CellIterator it ( *this );
	it.data = (void*) new Mesh::data_for_iter_max_dim; // !
	it.re_set = Mesh::reset_iter_max_dim;
	it.ad_vance = Mesh::advance_iter_max_dim;
	it.va_lid = Mesh::valid_iter_max_dim;
	if ( this->is_positive() )
		it.dereference = Mesh::deref_iter_pos_max_dim;
	else
	{	assert ( this->hidden_reverse != NULL );
		it.dereference = Mesh::deref_iter_neg_max_dim;    }
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Cells &, const tag::OfDim &, const short int d, const tag::NotOriented & )

// iterator over Cells belonging to 'this' mesh
// code example :  auto it = msh.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	CellIterator it ( *this );
	Mesh::data_for_iter_given_dim *dat = new Mesh::data_for_iter_given_dim; // !
	dat->dim = d;
	it.data = (void*) dat;
  it.re_set = Mesh::reset_iter_given_dim;
  it.ad_vance = Mesh::advance_iter_given_dim;
	it.va_lid = Mesh::valid_iter_given_dim;
	it.dereference = Mesh::deref_iter_given_dim;
	return it;                                                }


inline CellIterator Mesh::iter_over ( const tag::Segments &, const tag::Along & )

// iterator over segments of 'this' mesh
// beware : 'this' mesh should be open (a polygonal line)
// code example :  auto it = msh.iter_over ( tag::segments, tag::along );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_seg;
	it.ad_vance = Mesh::advance_iter_contour;
	it.va_lid = Mesh::valid_iter_segs_along;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                            }


inline CellIterator Mesh::iter_over ( const tag::Segments &, const tag::Around & )

// iterator over segments of 'this' mesh
// beware : 'this' mesh should be a closed loop (a polygon)
// code example :  auto it = msh.iter_over ( tag::segments, tag::around );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around;
	it.ad_vance = Mesh::advance_iter_contour;
	it.va_lid = Mesh::valid_iter_around;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                  }


inline CellIterator Mesh::iter_over ( const tag::Vertices &, const tag::Along & )

// iterator over vertices of 'this' mesh
// beware : 'this' mesh should be open (a polygonal line)
// code example :  auto it = msh.iter_over ( tag::vertices, tag::along );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_ver;
	it.ad_vance = Mesh::advance_iter_contour;
	it.va_lid = Mesh::valid_iter_vertices_along;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                           }


inline CellIterator Mesh::iter_over ( const tag::Vertices &, const tag::Around & )

// iterator over vertices of 'this' mesh
// beware : 'this' mesh should be a closed loop (a polygon)
// code example :  auto it = msh.iter_over ( tag::vertices, tag::around );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around;
	it.ad_vance = Mesh::advance_iter_contour;
	it.va_lid = Mesh::valid_iter_around;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Segments &, const tag::Along &, const tag::Reverse & )

// iterator over segments of to 'this' mesh
// beware : 'this' mesh should be open (a polygonal line)
// code example :  auto it = msh.iter_over ( tag::segments, tag::along, tag::reverse );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_seg_rev;
	it.ad_vance = Mesh::advance_iter_contour_rev;
	it.va_lid = Mesh::valid_iter_segs_along;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                          }


inline CellIterator Mesh::iter_over
	( const tag::Segments &, const tag::Around &, const tag::Reverse & )

// iterator over segments of 'this' mesh
// beware : 'this' mesh should be a closed loop (a polygon)
// code example :  auto it = msh.iter_over ( tag::segments, tag::around, tag::reverse );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_rev;
	it.ad_vance = Mesh::advance_iter_contour_rev;
	it.va_lid = Mesh::valid_iter_around;
	it.dereference = Mesh::deref_iter_segs_contour;
	return it;                                                            }


inline CellIterator Mesh::iter_over
	( const tag::Vertices &, const tag::Along &, const tag::Reverse & )

// iterator over vertices of 'this' mesh
// beware : 'this' mesh should be open (a polygonal line)
// code example :  auto it = msh.iter_over ( tag::vertices, tag::along, tag::reverse );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_along_ver_rev;
	it.ad_vance = Mesh::advance_iter_contour_rev;
	it.va_lid = Mesh::valid_iter_vertices_along;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                            }


inline CellIterator Mesh::iter_over
	( const tag::Vertices &, const tag::Around &, const tag::Reverse & )

// iterator over vertices of 'this' mesh
// beware : 'this' mesh should be a closed loop (a polygon)
// code example :  auto it = msh.iter_over ( tag::vertices, tag::around, tag::reverse );
//                 for ( it.reset(); it.in_range(); it++ ) { do something with *it; }

{	assert ( this->dim == 1 );
	CellIterator it ( *this );
	Mesh::data_for_iter_contour *d = new Mesh::data_for_iter_contour; // !
	d->oriented_base = this;
	it.data = (void*) d;
	it.re_set = Mesh::reset_iter_around_rev;
	it.ad_vance = Mesh::advance_iter_contour_rev;
	it.va_lid = Mesh::valid_iter_around;
	it.dereference = Mesh::deref_iter_vertices_contour;
	return it;                                                           }


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


// Fields implement the idea of coordinates (position of points in space).
// They can also hold other numeric values like the solution of a PDE,
// or material coefficients.

// Classes to be used in programs are OneDimField and MultiDimFields.
// Less often, one may need to use the base class NumericField.
// Classes which the end user will probably never need are :
// BlockFieldBase, BlockField, ComposedField.

// OneDimFields and BlockFields share the same internal organization of data
// (pointers into heaps attached to Cells) but are syntactically different
// because the operator() returns a double in the former case and a PointerVector
// of doubles in the latter case.
// Class PointerVector is used to extract the coordinates of a cell from a field;
// the user can change the values of the field by changing this PointerVector.

// BlockFields and CompositeFields are internally very different but exhibit
// the same behaviour. MultiDimFields' purpose is to provide a unified interface
// (a wrapper) to BlockFields and CompositeFields.

// most methods are defined in 'Field.cpp'


class ManiFEM::NumericField

// just a container for BlockFieldBase and ComposedField

{	public :

	std::string type;
	// "one-dim" or "multi-dim" or "block" or "composed"
  
	std::list < RelationFuncField * > relations;

	virtual size_t size () const; // assert 0

	// static std::map <short int, size_t> unused_index;

	// inline NumericField ()  {  assert (0); };
	// although syntactically this is not an abstract class,
	// I wanted to forbid (at run time) instantiation of a NumericField object 
	// but this doesn't work because C++ calls ALLWAYS the constructor of the base class :-(
  
	// factory functions :

	// Builds a one-dimensional field (one component), living on cells of a certain type;
	// 'type' can be "points", "segments", "2d cells" or "3d cells".
	// example of use :  auto &x = NumericField::one_dim ("lives on", "points");
	static OneDimField & one_dim ( std::string, std::string type );

	// Builds a field with n components, living on cells of a certain type;
	// 'type' can be "points", "segments", "2d cells" or "3d cells".
	// example of use :  auto &xy = NumericField::multi_dim ("dimension", 2, "lives on", "points");
	static MultiDimField & multi_dim ( std::string, size_t n, std::string, std::string type );

	// juxtaposition of two fields
	virtual MultiDimField & operator&& ( NumericField &f );
	virtual MultiDimField & join_when_first_is_block ( BlockFieldBase *first );
	virtual MultiDimField & join_when_first_is_composed ( ComposedField *first );
	virtual MultiDimField & join_when_first_is_multi ( MultiDimField *first );

	// projection (get one coordinate from a pair or a triplet)
	virtual OneDimField & operator[] (size_t i);
	
	// interpolation on the Mesh
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double, Cell&, double, Cell&, double );
  
}; // end of class ManiFEM::NumericField

//////////////////////////////////////////////////////////////////////////////

class ManiFEM::BlockFieldBase : public NumericField

// it's like a map from cells to (vectors of) real numbers
// One-dimensional Fields are specialized in OneDimField,
// while Fields with several components are specialized in BlockField.

{	public :

	// a Field is supposed to live on cells of a certain dimension
	short int lives_on_cells_of_dim;
	
	// within a cell's real_heap, 'this' field will occupy a certain set of indices
	size_t index_min, index_max_plus_one;

	std::set < NumericField* > is_part_of;
	// handled correctly by operator[]
	// not yet by operator&&

	// inline BlockFieldBase ()  {  assert (0); };
	// although syntactically this is not an abstract class,
	// I wanted to forbid (at run time) instantiation of a BlockFieldBase object
	// but this doesn't work because C++ calls ALLWAYS the constructor of the base class :-(
  
	virtual size_t size () const; // assert 0

	// coordinates of a cell
	PointerVector<double> & operator() ( Cell& );

	// projection (get one coordinate from a pair or a triplet)
	virtual OneDimField & operator[] ( size_t i ) = 0;
	
	// operator && is used to "concatenate" two fields
	virtual MultiDimField & operator&& ( NumericField &second );
	virtual MultiDimField & join_when_first_is_block ( BlockFieldBase *first );
	virtual MultiDimField & join_when_first_is_composed ( ComposedField *first );
	virtual MultiDimField & join_when_first_is_multi ( MultiDimField *first );

	// interpolation on the Mesh
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double, Cell&, double, Cell&, double );

}; // end of class ManiFEM::BlockFieldBase

//////////////////////////////////////////////////////////////////////////////

class ManiFEM::OneDimField : public BlockFieldBase

// it's like a map from cells to real numbers
// OneDimFields are built using the factory function NumericField::one_dim,
// or through the projection operator []
// This class is almost identical to BlockField; the main difference is
// the operator () which returns a real number (a double).

{	public :

	// constructors :
  
	inline OneDimField ()  {  type = "one-dim";  };
  
	inline OneDimField ( size_t dim_cells ) : OneDimField ()
	// consider using the factory function NumericField::one_dim instead,
	// except if you work with cells of dimension higher than 3
	{	this->lives_on_cells_of_dim = dim_cells;
		size_t & ind = Cell::real_heap_size [dim_cells];
		this->index_min = ind;
		ind ++;
		this->index_max_plus_one = ind;                         }

	virtual size_t size() const;

	// coordinate of a cell
	inline double & operator() ( Cell & cll ) const
	{	assert ( lives_on_cells_of_dim == cll.dim );
		assert ( cll.real_heap.size() > index_min );
		return cll.real_heap [index_min];              }

	// abstract coordinate on an arbitrary cell
	inline FieldOnAnyCell & operator() ( AnyCell & );
	// defined after class ManiFEM::FieldOnAnyCell

	// projection (get one coordinate from a pair or a triplet)
	virtual OneDimField & operator[] ( size_t i );

}; // end of class ManiFEM::OneDimField

//////////////////////////////////////////////////////////////////////////////

class ManiFEM::BlockField : public BlockFieldBase

// it's like a map from cells to vectors of real numbers
// BlockFields are built using the factory function NumericField::multi_dim,
// This class is almost identical to OneDimField; the main difference is
// the operator () which returns a vector <double>.

{	public :

	// constructors :
	// do not use directly, use instead factory functions
	// NumericField::one_dim and NumericField::multi_dim
  
	inline BlockField ()  {  type = "block";  };
  
	inline BlockField ( size_t n, short int dim_cells )
	{	// assert ( n > 1 );
		this->lives_on_cells_of_dim = dim_cells;
		size_t & ind = Cell::real_heap_size [dim_cells];
		this->index_min = ind;
		ind += n;
		this->index_max_plus_one = ind;                   }
	
	virtual size_t size() const;

	// projection (get one coordinate from a pair or a triplet)
	virtual OneDimField & operator[] ( size_t i );
};

//////////////////////////////////////////////////////////////////////////////

class ManiFEM::ComposedField : public NumericField

// A ComposedField is just a set of Fields (pair, triplet, etc)
// Used mainly to pack together several Fields through the operator &&.
// What is the difference between a ComposedField and a BlockField ?
// A ComposedField could be a pair of fields for instance.
// A corresponding BlockField would be a field of pairs of real numbers.
// BlockFields should be faster, but ComposedFields are often handy.
// BlockFields are built using the factory function BlockFieldBase::multi_dim,
// while ComposedFields are obtained through the operator &&.
// On both ComposedFields and BlockFields, the operator () returns a vector of real numbers.

// There is a special case when a BlockField cannot be used instead of a ComposedField :
// when we want to pack together two (or more) fields living on cells of different dimensions.

{	public :
	
	std::vector <BlockFieldBase*> fields;

	// constructor
	// do not use directly, use instead factory functions
	// NumericField::one_dim and NumericField::multi_dim
	inline ComposedField  (size_t n) : fields(n)
	{	type = "composed";  };

	virtual size_t size() const;

	// coordinates of a cell
	PointerVector<double> & operator() (Cell&);
	
	// projection (get one coordinate from a pair or a triplet)
	virtual OneDimField & operator[] (size_t i);
	
	// operator && is used to "concatenate" two fields
	virtual MultiDimField & operator&& (NumericField &f);
	virtual MultiDimField & join_when_first_is_block (BlockFieldBase *first);
	virtual MultiDimField & join_when_first_is_composed (ComposedField *first);
	virtual MultiDimField & join_when_first_is_multi (MultiDimField *first);

	// interpolation on the Mesh
	virtual void interpolate (Cell&, Cell&, double, Cell&, double);
	virtual void interpolate (Cell&, Cell&, double, Cell&, double, Cell&, double);
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double, Cell&, double, Cell&, double );

}; // end of class ManiFEM::ComposedField

//////////////////////////////////////////////////////////////////////////////

class ManiFEM::MultiDimField : public NumericField

// a unified interface for BlockField and ComposedField

{	public :

	void *base;

	virtual size_t size () const;

	// coordinates of a cell
	PointerVector<double> & (*restr) ( void*, Cell& );

	// projection (get one coordinate from a pair or a triplet)
	OneDimField & (*proj) (void*, size_t);

	// operator && is used to "concatenate" two fields
	MultiDimField & (*join_this_block) ( MultiDimField*, BlockFieldBase* );
	MultiDimField & (*join_this_composed) ( MultiDimField*, ComposedField* );
	MultiDimField & (*join_this_multi) ( MultiDimField*, MultiDimField* );
	MultiDimField & (*join_block_this) ( BlockFieldBase*, MultiDimField* );
	MultiDimField & (*join_composed_this) ( ComposedField*, MultiDimField* );
	MultiDimField & (*join_multi_this) ( MultiDimField*, MultiDimField* );
	
	// interpolation on the Mesh
	void (*interpolate2) ( MultiDimField*, Cell&, Cell&, double, Cell&, double );
	void (*interpolate3) ( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double );
	void (*interpolate4)
		( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double, Cell&, double );
	void (*interpolate6) ( MultiDimField*, Cell&, Cell&, double, Cell&, double,
     Cell&, double, Cell&, double, Cell&, double, Cell&, double );

	// (we are still in class MultiDimField) constructors :

	inline MultiDimField ()  { }
  
	inline MultiDimField ( size_t n, short int dim_cells ) : MultiDimField ()
	// consider using the factory function NumericField::multi_dim instead,
	// except if you work with cells of dimension higher than 3
 	{	assert ( n > 1 );
		type = "wrapped block";
		BlockField *b = new BlockField ( n, dim_cells );
		hidden_build (b);                                 }

	inline MultiDimField ( BlockFieldBase &b ) : MultiDimField (&b) { }
  
	inline MultiDimField ( BlockFieldBase *b ) : MultiDimField ()
	{	type = "wrapped block";
		hidden_build (b);        }
  
	inline MultiDimField ( ComposedField &b ) : MultiDimField (&b) { }
  
	inline MultiDimField ( ComposedField *b ) : MultiDimField ()
	{	type = "composed block";
		hidden_build (b);         }

	// used in constructors :

	void hidden_build ( BlockFieldBase* );
	void hidden_build ( ComposedField* );
	
	// (we are still in class MultiDimField) coordinates of a cell :

	static inline PointerVector<double> & restr_block ( void *b, Cell & cll )
	{	BlockFieldBase *base = (BlockFieldBase*) b;
		return (*base)(cll);                         }

	static inline PointerVector<double> & restr_composed ( void *b, Cell & cll )
	{	ComposedField *base = (ComposedField*) b;
		return (*base)(cll);                       }

	// projection (get one component from a pair or a triplet)

	static inline OneDimField & proj_block( void *b, size_t i )
	{	BlockField *base = (BlockField*) b;
		return (*base)[i];                       }

	static inline OneDimField & proj_composed ( void *b, size_t i )
	{	ComposedField *base = (ComposedField*) b;
		return (*base)[i];                             }

	// coordinates of a cell
	inline PointerVector<double> & operator() ( Cell & cll ) const
	{	return restr (base, cll);  }

	// projection (get one component from a pair or a triplet)
	virtual OneDimField & operator[] ( size_t i );

	// we are still in class MultiDimField
	// operator && is used to "concatenate" two fields
  
	virtual MultiDimField & operator&& ( NumericField &f );
  
	virtual MultiDimField & join_when_first_is_block ( BlockFieldBase *first );
	virtual MultiDimField & join_when_first_is_composed ( ComposedField *first );
	virtual MultiDimField & join_when_first_is_multi ( MultiDimField *first );
  
	static inline MultiDimField & join_wrappedblock_block
		( MultiDimField *first, BlockFieldBase *second  )
	{	return second->join_when_first_is_block ( (BlockFieldBase*) (first->base) );  }
  
	static inline MultiDimField & join_wrappedblock_composed
		( MultiDimField *first, ComposedField *second )
	{	return second->join_when_first_is_block ( (BlockFieldBase*) (first->base) );  }
	
	static inline MultiDimField & join_wrappedblock_multi
		( MultiDimField *first, MultiDimField *second )
	{	return second->join_when_first_is_block ( (BlockFieldBase*) (first->base) );  }

	static inline MultiDimField & join_wrappedcomposed_block
		( MultiDimField *first, BlockFieldBase *second )
	{	return second->join_when_first_is_composed ( (ComposedField*) (first->base) );  }
	
	static inline  MultiDimField & join_wrappedcomposed_composed
		( MultiDimField *first, ComposedField *second )
	{	return second->join_when_first_is_composed ( (ComposedField*) (first->base) );  }

	static inline MultiDimField & join_wrappedcomposed_multi
		( MultiDimField *first, MultiDimField *second )
	{	return second->join_when_first_is_composed ( (ComposedField*) (first->base) );  }

	static inline MultiDimField & join_block_wrappedblock
		( BlockFieldBase *first, MultiDimField *second )
	{	return ( (BlockFieldBase*) (second->base) ) -> join_when_first_is_block (first);  }

	static inline MultiDimField & join_composed_wrappedblock
		( ComposedField *first, MultiDimField *second )
	{	return ( (BlockFieldBase*) (second->base) ) -> join_when_first_is_composed (first);  }

	static inline MultiDimField & join_multi_wrappedblock
		( MultiDimField *first, MultiDimField *second )
	{	return ( (BlockFieldBase*) (second->base) ) -> join_when_first_is_multi (first);  }


	static MultiDimField & join_block_wrappedcomposed
		( BlockFieldBase *first, MultiDimField *second )
	{	return ( (ComposedField*) (second->base) ) -> join_when_first_is_block (first);  }

	static MultiDimField & join_composed_wrappedcomposed
		( ComposedField *first, MultiDimField *second )
	{	return ( (ComposedField*) (second->base) ) -> join_when_first_is_composed (first);  }

	static MultiDimField & join_multi_wrappedcomposed
		( MultiDimField *first, MultiDimField *second )
	{	return ( (ComposedField*) (second->base) ) -> join_when_first_is_multi (first);  }

	// we are still in class MultiDimField
	// interpolation on the Mesh

	virtual void interpolate ( Cell&, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double );
	virtual void interpolate ( Cell&, Cell&, double, Cell&, double, Cell&, double,
                                    Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_block_2
		( MultiDimField*, Cell&, Cell&, double, Cell&, double );
	static inline void interpolate_block_3
		( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_block_4
		( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_block_6 ( MultiDimField*, Cell&,
		Cell&, double, Cell&, double, Cell&, double,
	  Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_composed_2
		( MultiDimField*, Cell&, Cell&, double, Cell&, double );
	static inline void interpolate_composed_3
		( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_composed_4
		( MultiDimField*, Cell&, Cell&, double, Cell&, double, Cell&, double, Cell&, double );
	static inline void interpolate_composed_6 ( MultiDimField*, Cell&,
		Cell&, double, Cell&, double, Cell&, double,
	  Cell&, double, Cell&, double, Cell&, double );

}; // end of class MultiDimField

void inline MultiDimField::interpolate_block_2
(	MultiDimField *that, Cell & final, Cell & cell_1,
	double frac_1, Cell & cell_2, double frac_2        )
{	BlockFieldBase *b = (BlockFieldBase*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2 );  }

void inline MultiDimField::interpolate_block_3
(	MultiDimField *that, Cell & final, Cell & cell_1, double frac_1,
	Cell & cell_2, double frac_2, Cell & cell_3, double frac_3         )
{	BlockFieldBase *b = (BlockFieldBase*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3 );  }
	
void inline MultiDimField::interpolate_block_4
(	MultiDimField *that, Cell & final, Cell &cell_1, double frac_1, Cell & cell_2,
	double frac_2, Cell & cell_3, double frac_3, Cell & cell_4, double frac_4        )
{	BlockFieldBase *b = (BlockFieldBase*) that->base;
	b->interpolate
		( final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3, cell_4, frac_4 );  }
	
void inline MultiDimField::interpolate_block_6
(	MultiDimField *that, Cell & final, Cell &cell_1, double frac_1, Cell & cell_2,
	double frac_2, Cell & cell_3, double frac_3, Cell & cell_4, double frac_4,
	Cell & cell_5, double frac_5, Cell & cell_6, double frac_6                   )
{	BlockFieldBase *b = (BlockFieldBase*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2,
		cell_3, frac_3, cell_4, frac_4, cell_5, frac_5, cell_6, frac_6 );  }
	
void inline MultiDimField::interpolate_composed_2
(	MultiDimField *that, Cell & final, Cell & cell_1, double frac_1,
	Cell & cell_2, double frac_2                                     )
{	ComposedField *b = (ComposedField*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2 );  }
	
void inline MultiDimField::interpolate_composed_3
(	MultiDimField *that, Cell & final, Cell & cell_1, double frac_1, Cell &cell_2,
	double frac_2, Cell & cell_3, double frac_3                                     )
{	ComposedField *b = (ComposedField*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3 );  }
	
void inline MultiDimField::interpolate_composed_4
(	MultiDimField *that, Cell & final, Cell & cell_1, double frac_1, Cell &cell_2,
	double frac_2, Cell & cell_3, double frac_3, Cell & cell_4, double frac_4        )
{	ComposedField *b = (ComposedField*) that->base;
	b->interpolate
		( final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3, cell_4, frac_4 );  }

void inline MultiDimField::interpolate_composed_6
(	MultiDimField *that, Cell & final, Cell &cell_1, double frac_1, Cell & cell_2,
	double frac_2, Cell & cell_3, double frac_3, Cell & cell_4, double frac_4,
	Cell & cell_5, double frac_5, Cell & cell_6, double frac_6                   )
{	ComposedField *b = (ComposedField*) that->base;
	b->interpolate ( final, cell_1, frac_1, cell_2, frac_2,
		cell_3, frac_3, cell_4, frac_4, cell_5, frac_5, cell_6, frac_6 );  }
	

////////////////////////////////////////////////////////////////////////////////////

// Class PointerVector is used to extract the coordinates of a cell from a field;
// the user can change the values of the field by changing this PointerVector.
// It keeps an array of pointers towards the locations of the values.

template <typename VT> class ManiFEM::PointerVector

{	public :

	size_t size;

	// constructor
	inline PointerVector (size_t n) : pointers(n)
	{	size = n;  }

	// 'pointers' is an array of pointers towards the locations of the values.
	std::vector <VT*> pointers;

	inline VT & operator[] (size_t i) const
	{	return *(pointers[i]);  }

	// conversion operator used in manifolds (inner_prod, distance, projection, rotation)
	// not very efficient; is there another solution ?
	using vector_VT = std::vector<VT>;
	operator vector_VT () const
	{	vector_VT res (this->size);
		for ( size_t i = 0; i < this->size; i++ )
			res[i] = (*this)[i];
		return res;                                }

};


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

class ManiFEM::FieldOnAnyCell

// A NumericField applied to AnyCell
// used in abstract assertions of the type
// "on every point, this specific relation happens"

{	public:

	AnyCell * cell;
	OneDimField * field;

	inline FieldOnAnyCell ( AnyCell & ac )
	{	cell = & ac;  };
	
}; // end of class ManiFEM::FieldOnAnyCell

inline FieldOnAnyCell & OneDimField::operator() ( AnyCell & ac )
{	FieldOnAnyCell * f_ac = new FieldOnAnyCell (ac);
	f_ac->field = this;
	return *f_ac;                                       }

////////////////////////////////////////////////////////////////////////////////

class ManiFEM::FunctionOnAnyCell

// an FE_variable applied to AnyCell
// used in abstract assertions of the type
// "on every point, this specific relation happens"

{	public:

	AnyCell * cell;
	FunctionOnMesh::baseFunction * func;

	inline FunctionOnAnyCell ( FunctionOnMesh::baseFunction * f, AnyCell * ac )
	{	func = f;  cell = ac;  };
	
	inline void operator== ( FieldOnAnyCell & f_ac );
	inline void operator== ( double d );
	// defined after class FunctionOnMesh::Function
};

////////////////////////////////////////////////////////////////////////////////

class ManiFEM::RelationFunc

{	public:

	FunctionOnMesh::baseFunction * func;
	AnyCell * cell;
	std::string type;

	virtual std::pair<bool,double> operator() ( Cell & ) = 0;

};

////////////////////////////////////////////////////////////////////////////////

class ManiFEM::RelationFuncField : public RelationFunc

// class used in abstract assertions of the type
// "on every point, this specific relation happens"
// models a relation between an FunctionOnMesh::Function and a NumericField (OneDimField)

{	public:

	OneDimField * field;

	// constructor

	inline RelationFuncField ()
	{	type = "field";  }
	
	inline virtual std::pair<bool,double> operator() ( Cell & cll )
	{	std::pair <bool,double> p;
		p.first = this->field->lives_on_cells_of_dim == cll.dim;
		if ( p.first ) p.second = (*(this->field))(cll);
		return p;                                                      }
};

////////////////////////////////////////////////////////////////////////////////

class ManiFEM::RelationFuncConst : public RelationFunc

// class used in abstract assertions of the type
// "on every point, this specific relation happens"
// models a relation between two functions
// for instance, between an unknown and a given function
// or between a test function and a constant (typically zero)

{	public:

	double val;
	
	// constructor

	inline RelationFuncConst ()
	{	type = "const";  }
	
	inline virtual std::pair<bool,double> operator() ( Cell & )
	{	std::pair <bool,double> p ( 1, this->val );
		return p;                               }

};

////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::baseFunction

// just a container for classes like Constant, Sum, Product and many others
// but objects of 'baseFunction' type will also appear

{	public :

	std::string type, name;
	Manifold * manifold;
	
	std::list < RelationFunc * > relations;

	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );

	inline FunctionOnMesh::combinIntegrals & integrate ( Mesh & msh );
	// defined after class FunctionOnMesh::combinIntegrals

	inline double integrate ( Mesh & msh, Integrator & integ );
	// defined after class ManiFEM::Integrator

	inline void operator== ( double );
	// defined after class ManiFEM::FiniteElement

	virtual double operator() ( Cell & );

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	s.insert(this);  }
	
	inline virtual void seek_content (	VariationalProblem * vp ) { }

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{ return this;  }

	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	if ( this == x ) return y;
		return this;                 }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{ return this;  }

	virtual std::string repr ( std::string from = "" );
	
}; // end of class FunctionOnMesh::baseFunction

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Function

// just a a thin wrapper around FunctionOnMesh::baseFunction
// I guess this is what they call "delegation" ...
// useful for arithmetic operators and for vector functions

{	public :

	std::vector < FunctionOnMesh::baseFunction * > components;

	// The domain of a function could be a Mesh or a Manifold.
	void * domain;
	std::string domain_type;

	// constructors
	inline Function ( size_t d = 1 ) : components(d) {}
	inline Function ( FunctionOnMesh::baseFunction & f ) : Function (1)
	{	components[0] = &f;  }
	inline Function ( std::vector < FunctionOnMesh::baseFunction * > comp )
		: components ( comp ) { }

	inline FunctionOnMesh::Function & operator[] ( size_t i )
	{	return * ( new FunctionOnMesh::Function ( *(components[i]) ) );  }

	inline FunctionOnMesh::Function & operator+ ( double b ) const
	{	return operator+ ( FunctionOnMesh::constant(b) );  }

	inline FunctionOnMesh::Function & operator+
		( const FunctionOnMesh::Function & b ) const
	{	assert ( components.size() == 1 );
		assert ( b.components.size() == 1 );
		FunctionOnMesh::Function * new_fe = new FunctionOnMesh::Function;
		new_fe->components[0] =
			FunctionOnMesh::sum ( components[0], b.components[0] );
		return *new_fe;                                                     }

	inline FunctionOnMesh::Function & operator- ( ) const // unary
	{	assert ( components.size() == 1 );
		FunctionOnMesh::Function * new_fe = new FunctionOnMesh::Function;
		new_fe->components[0] = FunctionOnMesh::product
			( FunctionOnMesh::minus_one.components[0], components[0] ) ;
		return *new_fe;                                                    }

	inline FunctionOnMesh::Function & operator- ( double b ) const
	{	return operator+ ( -b );  }

	inline FunctionOnMesh::Function & operator-
		( const FunctionOnMesh::Function & b ) const
	{ return operator+ ( -b ); }

	inline FunctionOnMesh::Function & operator* ( double b ) const
	{	return operator* ( FunctionOnMesh::constant(b) );  }

	inline FunctionOnMesh::Function & operator*
		( const FunctionOnMesh::Function & b ) const
	{	assert ( components.size() == 1 );
		assert ( b.components.size() == 1 );
		FunctionOnMesh::Function * new_fe = new FunctionOnMesh::Function;
		new_fe->components[0] = FunctionOnMesh::product ( components[0], b.components[0] );
		return *new_fe;                                                                       }

	inline FunctionOnMesh::Function & operator/ ( double b ) const
	{	return operator/ ( FunctionOnMesh::constant(b) );  }

	inline FunctionOnMesh::Function & operator/
		( const FunctionOnMesh::Function & b ) const
	{	assert ( components.size() == 1 );
		assert ( b.components.size() == 1 );
		FunctionOnMesh::Function * new_fe = new FunctionOnMesh::Function;
		new_fe->components[0] = FunctionOnMesh::fraction ( components[0], b.components[0] );
		return *new_fe;                                                                        }

	inline FunctionOnMesh::Function & deriv
		( const FunctionOnMesh::Function & var )
	// derivative with respect to 'var'
	{	assert ( components.size() == 1 );
		assert ( var.components.size() == 1 );
		FunctionOnMesh::Function * new_fe = new FunctionOnMesh::Function;
		new_fe->components[0] =  components[0]->deriv(var.components[0]);
		return * new_fe;                                                    }

	// we are still in class FunctionOnMesh::Function

	inline FunctionOnMesh::Function & operator&&
		( const FunctionOnMesh::Function & b ) const
	{	FunctionOnMesh::Function * res = new FunctionOnMesh::Function (0);
		size_t n = components.size() + b.components.size();
		res->components.reserve ( n );
		res->components.insert
			( res->components.end(), components.begin(), components.end() );
		res->components.insert
			( res->components.end(), b.components.begin(), b.components.end() );
		return * res;                                                            }
	
	inline void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	std::vector<FunctionOnMesh::baseFunction*>::iterator it;
		for ( it = components.begin(); it != components.end(); it++ )
			(*it)->seek_vars (s);                                        }

	inline FunctionOnMesh::Function & compose_with_inverse_of
		( FunctionOnMesh::Mapping & T );
	// defined after class FunctionOnMesh::Composed
	
	inline FunctionOnMesh::combinIntegrals & integrate ( Mesh & msh );
	// defined after class FunctionOnMesh::combinIntegrals

	inline double integrate ( Mesh & msh, Integrator & integ );
	// defined after class ManiFEM::Integrator

	inline void operator== ( FunctionOnMesh::Function & );
	inline void operator== ( double );
	// defined after class ManiFEM::FiniteElement

	double operator() ( Cell& );

	inline FunctionOnAnyCell & operator() ( AnyCell & ac )
	{	assert ( components.size() == 1 );
		FunctionOnAnyCell * f_ac =
			new FunctionOnAnyCell ( components[0], &ac );
		return *f_ac;                                       }

	virtual std::string repr ( std::string from = "" );
	
}; // end of class FunctionOnMesh::Function

// non-member operators :

inline FunctionOnMesh::Function & operator+
	( double c, const FunctionOnMesh::Function & b )
{ return FunctionOnMesh::constant(c) + b; }

inline FunctionOnMesh::Function & operator-
	( double c, const FunctionOnMesh::Function & b )
{	return FunctionOnMesh::constant(c) - b; }

inline FunctionOnMesh::Function & operator*
	( double c, const FunctionOnMesh::Function & b )
{ return FunctionOnMesh::constant(c) * b; }

inline FunctionOnMesh::Function & operator/
	( double c, const FunctionOnMesh::Function & b )
{ return FunctionOnMesh::constant(c) / b; }

// other operators :

inline void FunctionOnAnyCell::operator== ( FieldOnAnyCell & f_ac )
{	// keep record of an equality between 'this' and f_ac
	RelationFuncField * rel = new RelationFuncField;
	rel->func = this->func;
	rel->field = f_ac.field;
	assert ( f_ac.cell == this->cell );
	rel->cell = f_ac.cell;
	this->func->relations.push_back (rel);
	f_ac.field->relations.push_back (rel);             }

inline void FunctionOnAnyCell::operator== ( double d )
{	// keep record of an equality between 'this' and d
	RelationFuncConst * rel = new RelationFuncConst;
	rel->func = this->func;
	rel->cell = this->cell;
	rel->val  = d;
	this->func->relations.push_back (rel);               }

// other stuff

inline FunctionOnMesh::baseFunction * FunctionOnMesh::baseFunction::deriv
	( FunctionOnMesh::baseFunction * var )
// the derivative with respect to itself should be one
// the derivative with respect to any other variable should be zero
{	// assert that the manifold is the same ... not easy !
	if ( this == var ) return FunctionOnMesh::one.components[0];
	return FunctionOnMesh::zero.components[0];                      }


//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Constant
	: public FunctionOnMesh::baseFunction

// a constant function

{	public :

	double c;

	// constructor (do not use; use factory functions instead)
	inline Constant (double cc) : FunctionOnMesh::baseFunction ()
	{	type = "constant"; c = cc;  }

	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var )
	// derivative with respect to 'var'
	{	return FunctionOnMesh::zero.components[0];  }

	inline virtual double operator() ( Cell& )
	{	return c;  }

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & )	{ }
	
	inline virtual void seek_content (	VariationalProblem * vp ) { }

	// 'replace', 'replace_var' and 'compute_derivatives'
	// inherited from class FunctionOnMesh::baseFunction

	virtual std::string repr ( std::string from = "" );

}; // end of class FunctionOnMesh::Constant

// factory function :

inline FunctionOnMesh::Function & FunctionOnMesh::constant (double cc)
{	if ( cc == 0. ) return FunctionOnMesh::zero;
	if ( cc == 1. ) return FunctionOnMesh::one;
	if ( cc == -1. ) return FunctionOnMesh::minus_one;
	FunctionOnMesh::baseFunction * c =
		(FunctionOnMesh::baseFunction*) new FunctionOnMesh::Constant (cc);
	return * ( new FunctionOnMesh::Function (*c) );                       };

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Unknown
	: public FunctionOnMesh::baseFunction

// a function related to a NumericField, to be computed
// specifies how the function relates to the field
// for instance, the field keeps the values of the function at vertices
// between vertices the function varies linearly

{	public:

	// Manifold * manifold;
	NumericField * field;
	FunctionOnMesh::Test * test;
	
	// constructor
	inline Unknown ( NumericField & f ) :
		FunctionOnMesh::baseFunction ()
	{	size_t d = f.size();
		type = "unknown";
		field = & f;
		// manifold = new Manifold (d);
  }

	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );
	// derivative with respect to 'var'
	// defined after class FunctionOnMesh::derivative

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & )
	{	assert (0);  }
	
	inline virtual void seek_content ( VariationalProblem * vp );
	// defined after class ManiFEM::VariationalProblem

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	if ( f == this ) return g;
		return this;                }
	
	// 'replace_var' and 'repr' inherited from class FunctionOnMesh::baseFunction
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{ assert (0);  }

}; // end of class FunctionOnMesh::Unknown

// factory function

inline FunctionOnMesh::Function & FunctionOnMesh::unknown
	( NumericField & field, const std::string type )
{	assert ( field.size() == 1 );
	OneDimField & f = * ( (OneDimField*) &field );
	// we would love to build Functions from Fields with more than
	// one component but for the time being we don't know how ... :-(
	FunctionOnMesh::Unknown * u = new FunctionOnMesh::Unknown ( field );
	assert ( type == "Lagrange degree one");
	// we provide a NULL mesh, which means "for all cells, in any mesh" :
	AnyCell & cll = * ( new AnyCell ( (Mesh*) NULL, 0 ) );
	FunctionOnMesh::Function * uu = new FunctionOnMesh::Function (*u);
	(*uu) (cll) == f (cll);
	return *uu;                                                              }

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Test
	: public FunctionOnMesh::baseFunction

// a function to be used as test in variational formulations
// makes pair with an unknown function
// no field associated

{	public:

	FunctionOnMesh::Unknown * unknown;

	// constructor
	inline Test ( FunctionOnMesh::Unknown * u ) :
		FunctionOnMesh::baseFunction ()
	{	unknown = u;  u->test = this;   }

	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );
	// derivative with respect to 'var'
	// defined after class FunctionOnMesh::derivative

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & )
	{	assert (0);  }
	
	inline virtual void seek_content (	VariationalProblem * vp );
	// defined after class ManiFEM::VariationalProblem

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	if ( f == this ) return g;
		return this;                }
	
	// 'replace_var' and 'repr' inherited from class FunctionOnMesh::baseFunction

	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{ assert (0);  }

}; // end of class FunctionOnMesh::Test

// factory function

inline FunctionOnMesh::Function & FunctionOnMesh::test
	( FunctionOnMesh::Function & u )
{	assert ( u.components.size() == 1 );
	FunctionOnMesh::Test * t = new FunctionOnMesh::Test
		( (FunctionOnMesh::Unknown*) u.components[0] );
	return * ( new FunctionOnMesh::Function (*t) );         }

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Sum
	: public FunctionOnMesh::baseFunction

// sum between FunctionOnMesh::baseFunction

{	public :

	std::list < FunctionOnMesh::baseFunction * > terms;

	// constructors (do not use; use factory functions instead)
	
	inline Sum ( ) : FunctionOnMesh::baseFunction ()
	{	type = "sum";  }
		
	inline Sum ( std::list < FunctionOnMesh::baseFunction * > l ) :
		Sum ()  { terms = l;  }

	// derivative with respect to 'var'
	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );

	inline virtual double operator() ( Cell & cll )
	{	double res = 0.;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			res += (*it)->operator() (cll);
		return res;                                                         	}

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			(*it)->seek_vars (s);                                         }
	
	inline virtual void seek_content ( VariationalProblem * vp )
	{	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			(*it)->seek_content ( vp );                                   }

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Sum * res = new FunctionOnMesh::Sum;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			res->terms.push_back ( (*it)->replace ( f, g ) );
		return res;                                          }
	
	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	FunctionOnMesh::Sum * res = new FunctionOnMesh::Sum;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			res->terms.push_back ( (*it)->replace_var ( x, y ) );
		return res;                                          }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{	FunctionOnMesh::Sum * res = new FunctionOnMesh::Sum;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = terms.begin(); it != terms.end(); it++ )
			res->terms.push_back ( (*it)->compute_derivatives() );
		return res;                                          }
	
	virtual std::string repr ( std::string from = "" );
	
}; // end of class FunctionOnMesh::Sum

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Sum::deriv
	( FunctionOnMesh::baseFunction * var )
{	std::list < FunctionOnMesh::baseFunction * > l;
	for ( std::list<FunctionOnMesh::baseFunction*>::iterator it =
				terms.begin(); it != terms.end(); it++ )
	{	l.push_back ( (*it)->deriv(var) );    };
	return FunctionOnMesh::sum (l);                            }

// factory functions :

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::sum
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
{	std::list < FunctionOnMesh::baseFunction * > ter;
	ter.push_back(f);
	ter.push_back(g);
	return FunctionOnMesh::sum (ter);              }

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::sum
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g,
  FunctionOnMesh::baseFunction * h )
{	std::list < FunctionOnMesh::baseFunction * > ter;
	ter.push_back(f);
	ter.push_back(g);
	ter.push_back(h);
	return FunctionOnMesh::sum (ter);               }

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::sum
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g,
  FunctionOnMesh::baseFunction * h, FunctionOnMesh::baseFunction * j )
{	std::list < FunctionOnMesh::baseFunction * > ter;
	ter.push_back(f);
	ter.push_back(g);
	ter.push_back(h);
	ter.push_back(j);
	return FunctionOnMesh::sum (ter);               }

////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Product
	: public FunctionOnMesh::baseFunction

// product between FunctionOnMesh::baseFunction

{	public :

	std::list < FunctionOnMesh::baseFunction * > factors;

	// constructors (do not use; use factory functions instead)
	
	inline Product () : FunctionOnMesh::baseFunction ()
	{	type = "product";  }
	
	inline Product ( std::list < FunctionOnMesh::baseFunction * > l ) :
		Product ()  { factors = l;  }

	// derivative with respect to 'var'
	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );

	inline virtual double operator() ( Cell & cll )
	{	double res = 1.;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = factors.begin(); it != factors.end(); it++ )
			res *= (*it)->operator() (cll);
		return res;                                                         	}
	
	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = this->factors.begin(); it != this->factors.end(); it++ )
			(*it)->seek_vars (s);                                         }
	
	inline virtual void seek_content ( VariationalProblem * vp )
	{	for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = factors.begin(); it != factors.end(); it++ )
			(*it)->seek_content ( vp );                                      }

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Product * res = new FunctionOnMesh::Product;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = factors.begin(); it != factors.end(); it++ )
			res->factors.push_back ( (*it)->replace ( f, g ) );
		return res;                                                   }
	
	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	FunctionOnMesh::Product * res = new FunctionOnMesh::Product;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = factors.begin(); it != factors.end(); it++ )
			res->factors.push_back ( (*it)->replace_var ( x, y ) );
		return res;                                                   }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{	FunctionOnMesh::Product * res = new FunctionOnMesh::Product;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator
	        it = factors.begin(); it != factors.end(); it++ )
			res->factors.push_back ( (*it)->compute_derivatives() );
		return res;                                                   }
	
	virtual std::string repr ( std::string from = "" );

}; // end of class FunctionOnMesh::Product

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Product::deriv
	( FunctionOnMesh::baseFunction * var )
{	std::list < FunctionOnMesh::baseFunction * > l;
	for ( std::list<FunctionOnMesh::baseFunction*>::iterator it = this->factors.begin();
				it != factors.end(); it++ )
	{	std::list < FunctionOnMesh::baseFunction * > ll;
		for ( std::list<FunctionOnMesh::baseFunction*>::iterator itt = this->factors.begin();
	        itt != factors.end(); itt++ )
			if ( it == itt ) ll.push_back ( (*itt)->deriv(var) );
			else ll.push_back ( *itt );
		l.push_back ( FunctionOnMesh::product (ll) );                                }
	return FunctionOnMesh::sum (l);                                                   }

// factory functions :

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::product
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
{	std::list < FunctionOnMesh::baseFunction * > fact;
	fact.push_back(f);
	fact.push_back(g);
	return FunctionOnMesh::product (fact);          }

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::product
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g,
  FunctionOnMesh::baseFunction * h )
{	std::list < FunctionOnMesh::baseFunction * > fact;
	fact.push_back(f);
	fact.push_back(g);
	fact.push_back(h);
	return FunctionOnMesh::product (fact);          }

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::product
( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g,
  FunctionOnMesh::baseFunction * h, FunctionOnMesh::baseFunction * j )
{	std::list < FunctionOnMesh::baseFunction * > fact;
	fact.push_back(f);
	fact.push_back(g);
	fact.push_back(h);
	fact.push_back(j);
	return FunctionOnMesh::product (fact);          }

////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Fraction
	: public FunctionOnMesh::baseFunction

// product between FunctionOnMesh::baseFunction

{	public :

	FunctionOnMesh::baseFunction * up, * down;

	// constructors (do not use; use factory functions instead)

	inline Fraction () : FunctionOnMesh::baseFunction ()
	{	type = "fraction";  }
	
	inline Fraction ( FunctionOnMesh::baseFunction * f,
	                  FunctionOnMesh::baseFunction * g )
		: Fraction ()  { up = f; down = g; };

	// derivative with respect to 'var'
	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );
		
	inline virtual double operator() ( Cell & cll )
	{	return up->operator()(cll) /down->operator()(cll); 	}
	
	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	up->seek_vars (s);  down->seek_vars (s);  }
	
	inline virtual void seek_content ( VariationalProblem * vp )
	{	up->seek_content ( vp );  down->seek_content ( vp );  }

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Fraction * res = new FunctionOnMesh::Fraction;
		res->up   = up  ->replace ( f, g );
		res->down = down->replace ( f, g );
		return res;                                                      }
	
	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	FunctionOnMesh::Fraction * res = new FunctionOnMesh::Fraction;
		res->up   = up  ->replace_var ( x, y );
		res->down = down->replace_var ( x, y );
		return res;                                                      }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{	FunctionOnMesh::Fraction * res = new FunctionOnMesh::Fraction;
		res->up   = up  ->compute_derivatives();
		res->down = down->compute_derivatives();
		return res;                                                      }
	
	virtual std::string repr ( std::string from = "" );

}; // end of class FunctionOnMesh::Fraction

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Fraction::deriv
	( FunctionOnMesh::baseFunction * var )
{	FunctionOnMesh::baseFunction * a = FunctionOnMesh::product ( up->deriv(var), down );
	FunctionOnMesh::baseFunction * b = FunctionOnMesh::product ( up, down->deriv(var) );
	FunctionOnMesh::baseFunction * c = FunctionOnMesh::product
		( FunctionOnMesh::minus_one.components[0], b );
	FunctionOnMesh::baseFunction * dif = FunctionOnMesh::sum ( a, c );
	FunctionOnMesh::baseFunction * sq = FunctionOnMesh::product ( down, down );
	return FunctionOnMesh::fraction ( dif, sq );                                          }
		
// factory function :

inline FunctionOnMesh::baseFunction *	FunctionOnMesh::fraction
	( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
{
	if ( g->type == "constant" )
	{
		FunctionOnMesh::Constant * gg = (FunctionOnMesh::Constant*) g;
		double c = gg->c;
		assert ( c != 0. );
		gg = new FunctionOnMesh::Constant ( 1./c );
		g = (FunctionOnMesh::baseFunction*) gg;
		return FunctionOnMesh::product ( f, g );                         }
	return new FunctionOnMesh::Fraction ( f, g );                         }

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Power
	: public FunctionOnMesh::baseFunction

// a baseFunction elevated to a (constant) power

{	public :
  
	FunctionOnMesh::baseFunction * base;
	double exponent;
  
	// constructors (do not use; use factory functions instead)

	inline Power () : FunctionOnMesh::baseFunction ()
	{	type = "power";  }

	inline Power ( FunctionOnMesh::baseFunction * b, double e )
		: Power ()
	{	base = b; exponent = e; };

	inline virtual double operator() ( Cell & cll )
	{	return pow ( base->operator()(cll), exponent ); 	}
	
	// derivative with respect to 'var'
	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var )
	{	FunctionOnMesh::baseFunction * exp  =
			(FunctionOnMesh::baseFunction*) new FunctionOnMesh::Constant ( exponent );
		FunctionOnMesh::baseFunction * p = FunctionOnMesh::power ( base, exponent-1. );
		return  ( product ( exp, p, base->deriv ( var ) ) );                             };

	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & s )
	{	base->seek_vars (s);  }
	
	inline virtual void seek_content ( VariationalProblem * vp )
	{	base->seek_content ( vp );  }

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Power * res = new FunctionOnMesh::Power;
		res->base = base->replace ( f, g );
		res->exponent = exponent;
		return res;                                               }
	
	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	FunctionOnMesh::Power * res = new FunctionOnMesh::Power;
		res->base = base->replace_var ( x, y );
		res->exponent = exponent;
		return res;                                               }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{	FunctionOnMesh::Power * res = new FunctionOnMesh::Power;
		res->base = base->compute_derivatives();
		res->exponent = exponent;
		return res;                                               }
	
	virtual std::string repr ( std::string from = "" );

}; // end of class FunctionOnMesh::Power

// factory function :

inline FunctionOnMesh::baseFunction * FunctionOnMesh::power
	( FunctionOnMesh::baseFunction * b, double e )
{	// we should give a special treatment to the case when 'b' is a power
	if ( e == 1. ) return b;
	if ( e == 0. ) return FunctionOnMesh::one.components[0];
	if ( b == FunctionOnMesh::zero.components[0] ) return FunctionOnMesh::zero.components[0];
	if ( b == FunctionOnMesh::one.components[0] ) return FunctionOnMesh::one.components[0];
	if ( b->type == "power" )
	{	FunctionOnMesh::Power * bb = (FunctionOnMesh::Power*) b;
		e *= bb->exponent;
		b = bb->base;                   }
	return (FunctionOnMesh::baseFunction*) new FunctionOnMesh::Power (b, e);  };

////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Mapping
	: public FunctionOnMesh::Function

// a diffeomeorphism, a change of coordinates

{	public :

	std::set < FunctionOnMesh::baseFunction* > variables;

	FunctionOnMesh::baseFunction * det;
	std::map < FunctionOnMesh::baseFunction *,
	      std::map < FunctionOnMesh::baseFunction *,
	            FunctionOnMesh::baseFunction * > > jacobian, inv_jacobian;
	// The method "build_jacobian" computes symbolically the partial derivatives
	// and stores them in "jacobian". If this matrix is square, it also computes
	// the inverse matrix (as symbolic functions) and stores them in "inv_jacobian".

	// constructor
	inline Mapping ( std::vector < FunctionOnMesh::baseFunction * > comp )
		: FunctionOnMesh::Function ( comp ) { } //	{ components = comp; }

	void build_jacobian ();
	// The method "build_jacobian" computes symbolically the partial derivatives
	// and stores them in "jacobian". If this matrix is square, it also computes
	// the inverse matrix (as symbolic functions) and stores them in "inv_jacobian".

}; // end of class FunctionOnMesh::Mapping

// factory function

inline FunctionOnMesh::Mapping & FunctionOnMesh::mapping
	( FunctionOnMesh::Function & f )

{	// we create a new object, a Mapping, with the same components as f
	FunctionOnMesh::Mapping * new_f =
		new FunctionOnMesh::Mapping ( f.components );
	// destroy f :
	delete ( &f );
	// build list of variables :
	new_f->seek_vars ( new_f->variables );
	// assert ( new_f->variables.size() == new_f->components.size() );
	// no need to size the jacobian matrix - maps adapt themselves
	new_f->build_jacobian();
	return * new_f;                                                       }

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Composition
	: public FunctionOnMesh::baseFunction

// a function composed with the INVERSE of a Mapping

{	public :

	FunctionOnMesh::Function * func;
	FunctionOnMesh::Mapping * inv_transf;
	
	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var );

	inline virtual double operator() ( Cell & cll )
	{	return func->operator()(cll); 	}
	
	inline virtual void seek_vars ( std::set < FunctionOnMesh::baseFunction* > & )
	{	assert (0);  }
	
	inline virtual void seek_content ( VariationalProblem * vp )
	{	assert (0);  }

	// 'replace' and 'replace_var' inherited from class FunctionOnMesh::baseFunction

	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{ assert (0);  }

	std::string repr ( std::string from = "" );

}; // end of class FunctionOnMesh::Composition

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Composition::deriv
	( FunctionOnMesh::baseFunction * var )

{	std::map < FunctionOnMesh::baseFunction *,
	      FunctionOnMesh::baseFunction * > &
		column_of_jacobian = this->inv_transf->inv_jacobian.at(var);
	FunctionOnMesh::baseFunction * res = FunctionOnMesh::zero.components[0];
	std::set<FunctionOnMesh::baseFunction*>::iterator it;
	for ( it  = this->inv_transf->variables.begin();
	      it != this->inv_transf->variables.end();   it++ )
	{	FunctionOnMesh::baseFunction * interm_var = *it;
		FunctionOnMesh::baseFunction * temp = FunctionOnMesh::product
			( this->func->components[0]->deriv (interm_var),
			  column_of_jacobian.at(interm_var)              );
		res = FunctionOnMesh::sum ( res, temp );                       }
	return res;                                                                }
	
inline FunctionOnMesh::Function &
FunctionOnMesh::Function::compose_with_inverse_of
	( FunctionOnMesh::Mapping & T )

// when we compose 'this' with the inverse of T :
// check for the variables of 'this' (on which 'this' depends)
// do the same thing for T and make sure that they are the same
// then, simply link with T; the values of T will be the new variables

{	std::set < FunctionOnMesh::baseFunction* > vars_of_this, vars_of_T;
	this->seek_vars ( vars_of_this );
	T.seek_vars ( vars_of_T );
	assert ( vars_of_this == vars_of_T );
	assert ( vars_of_T.size() == T.components.size() );
	assert ( this->components.size() == 1 );
	FunctionOnMesh::Composition * f = new FunctionOnMesh::Composition;
	f->func = this;
	f->inv_transf = &T;
	return * ( new FunctionOnMesh::Function
	              ( * (FunctionOnMesh::baseFunction*) f ) );              }

//////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Derivative
	: public FunctionOnMesh::baseFunction

// delayed derivative of some expression 'base'
// usually this expression involves an unknown or a test function
	
{	public:

	FunctionOnMesh::baseFunction * base, * var;

	inline virtual FunctionOnMesh::baseFunction * deriv
		( FunctionOnMesh::baseFunction * var )
	{	assert(0);  }
	
	inline virtual void seek_content ( VariationalProblem * vp );
	// defined after class ManiFEM::VariationalProblem

	inline virtual FunctionOnMesh::baseFunction * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Derivative * res = new FunctionOnMesh::Derivative;
		res->base = base->replace ( f, g );
		res->var = var;
		return res;                                                         }
	
	inline virtual FunctionOnMesh::baseFunction * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{
		FunctionOnMesh::Derivative * res = new FunctionOnMesh::Derivative;
		res->base = base->replace_var ( x, y );
		if ( var == x ) res->var = y;
		else res->var = var;
		return res;                                                         }
	
	inline virtual FunctionOnMesh::baseFunction * compute_derivatives ()
	{	return base->deriv(var);  }

	inline virtual std::string repr ( std::string from = "" )
	{	return "deriv(" + base->repr() + "," + var->repr() + ")";  }
	
}; // end of class FunctionOnMesh::Derivative

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Unknown::deriv
	( FunctionOnMesh::baseFunction * var )
// derivative with respect to 'var'
{	FunctionOnMesh::Derivative * res = new FunctionOnMesh::Derivative;
	res->base = this;  res->var = var;
	return (FunctionOnMesh::baseFunction*) res;                           }

inline FunctionOnMesh::baseFunction * FunctionOnMesh::Test::deriv
	( FunctionOnMesh::baseFunction * var )
// derivative with respect to 'var'
{	FunctionOnMesh::Derivative * res = new FunctionOnMesh::Derivative;
	res->base = this;  res->var = var;
	return (FunctionOnMesh::baseFunction*) res;                           }

//////////////////////////////////////////////////////////////////////////////////
////////////////         integrals         ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::Integral

// delayed integral of a function

{	public :
	
	Mesh * domain;

	FunctionOnMesh::baseFunction * integrand;

	// constructor
	// inline Integral ( FunctionOnMesh::baseFunction * in, Mesh * msh )
	// {	domain = msh; integrand = in; }

	inline void seek_content ( VariationalProblem * vp )
	{	integrand->seek_content ( vp ); }

	inline FunctionOnMesh::Integral * replace
		( FunctionOnMesh::baseFunction * f, FunctionOnMesh::baseFunction * g )
	{	FunctionOnMesh::Integral * res = new FunctionOnMesh::Integral;
		res->integrand = integrand->replace ( f, g );
		return res;                                                      }

	inline FunctionOnMesh::Integral * replace_var
		( FunctionOnMesh::Function & x, FunctionOnMesh::Function & y )
	{	assert( x.components.size() == y.components.size() );
		size_t n = x.components.size();
		FunctionOnMesh::Integral * res = this;
		for ( size_t i = 0; i < n; i++ )
			res = res->replace_var ( x.components[i], y.components[i] );
		return res;                                                      }

	// replace_var should be deeply reviewed
	
	inline FunctionOnMesh::Integral * replace_var
		( FunctionOnMesh::baseFunction * x, FunctionOnMesh::baseFunction * y )
	{	FunctionOnMesh::Integral * res = new FunctionOnMesh::Integral;
		res->domain = domain;
		res->integrand = integrand->replace_var ( x, y );
		return res;                                                      }

	inline FunctionOnMesh::Integral * compute_derivatives ()
	{	FunctionOnMesh::Integral * res = new FunctionOnMesh::Integral;
		res->domain = domain;
		res->integrand = integrand->compute_derivatives();
		return res;                                                     }

	inline virtual std::string repr ( std::string from = "" )
	{	return "int(" + integrand->repr() + ")";  }
	
}; // end of class FunctionOnMesh::Integral

////////////////////////////////////////////////////////////////////////////////

template <class U, class V> class ManiFEM::Equality
{	public :
	U * lhs; V * rhs;
	// constructor
	inline Equality () { }
	inline Equality ( U * u, V * v )
	{	lhs = u;  rhs = v;  }
};

////////////////////////////////////////////////////////////////////////////////

class FunctionOnMesh::combinIntegrals

// linear combination of integrals, not necessarily on the same domain

{	public :

	std::list < std::pair < double, FunctionOnMesh::Integral * > > terms;

	// constructors
	inline combinIntegrals () : terms() { };
	inline combinIntegrals ( FunctionOnMesh::baseFunction * f, Mesh * msh )
		: combinIntegrals ()
	{	terms.emplace_front();
		std::list < std::pair < double, FunctionOnMesh::Integral * > >
			::iterator lib = terms.begin();
		lib->first = 1.;
		lib->second = new FunctionOnMesh::Integral;
		lib->second->domain = msh;
		lib->second->integrand = f;                               };

	// arithmetic

	inline FunctionOnMesh::combinIntegrals &
		operator+ ( const FunctionOnMesh::combinIntegrals & b ) const;

	inline FunctionOnMesh::combinIntegrals &
		operator- ( const FunctionOnMesh::combinIntegrals & b ) const;

	inline FunctionOnMesh::combinIntegrals &
		operator- ( ); // unary

//	inline Equality < FunctionOnMesh::combinIntegrals,
//	                  FunctionOnMesh::combinIntegrals >
//		& operator== ( FunctionOnMesh::combinIntegrals & rhs )
//	{	Equality < FunctionOnMesh::combinIntegrals,
//		           FunctionOnMesh::combinIntegrals >
//			* res = new Equality < FunctionOnMesh::combinIntegrals,
//			                       FunctionOnMesh::combinIntegrals>;
//		res->lhs = this;  res->rhs = &rhs;
//		return * res;                                                }
	inline VariationalProblem
		& operator== ( FunctionOnMesh::combinIntegrals & rhs );
	// factory function for variational problems
	// defined after class ManiFEM::VariationalProblem

	inline virtual void seek_content ( VariationalProblem * vp )
	{	for ( std::list<std::pair<double,FunctionOnMesh::Integral*>>::iterator
	        it = this->terms.begin(); it != this->terms.end(); it++ )
			it->second->seek_content ( vp );                             }

}; // end of class FunctionOnMesh::combinIntegrals

inline FunctionOnMesh::combinIntegrals &
	FunctionOnMesh::combinIntegrals::operator+
	( const FunctionOnMesh::combinIntegrals & b ) const
{	FunctionOnMesh::combinIntegrals * res =
		new FunctionOnMesh::combinIntegrals ();
	for ( std::list < std::pair < double, FunctionOnMesh::Integral * > >
					::const_iterator it = terms.begin(); it != terms.end(); it++ )
		res->terms.push_back ( *it );
	for ( std::list < std::pair < double, FunctionOnMesh::Integral * > >
		::const_iterator it = b.terms.begin(); it != b.terms.end(); it++ )
		res->terms.push_back ( *it );
	return *res;                                                         }

inline FunctionOnMesh::combinIntegrals &
	FunctionOnMesh::combinIntegrals::operator-
	( const FunctionOnMesh::combinIntegrals & b ) const
{	FunctionOnMesh::combinIntegrals * res =
		new FunctionOnMesh::combinIntegrals ();
	for ( std::list < std::pair < double, FunctionOnMesh::Integral * > >
		::const_iterator it = terms.begin(); it != terms.end(); it++ )
		res->terms.push_back ( *it );
	for ( std::list < std::pair < double, FunctionOnMesh::Integral * > >
		::const_iterator it = b.terms.begin(); it != b.terms.end(); it++ )
	{	res->terms.push_front ( *it );  res->terms.begin()->first *= -1.;  }
	return *res;                                                          }

inline FunctionOnMesh::combinIntegrals &
	FunctionOnMesh::combinIntegrals::operator- ( ) // unary
{	for ( std::list < std::pair < double, FunctionOnMesh::Integral * > >
		::iterator it = terms.begin(); it != terms.end(); it++ )
		it->first *= -1.;
	return *this;                                                          }

inline FunctionOnMesh::combinIntegrals &
	FunctionOnMesh::baseFunction::integrate ( Mesh & msh )
{	return * ( new FunctionOnMesh::combinIntegrals ( this, &msh ) ); } 

inline FunctionOnMesh::combinIntegrals &
	FunctionOnMesh::Function::integrate ( Mesh & msh )
{	assert ( this->components.size() == 1 );
	return this->components[0]->integrate(msh); } 


/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
//      H       H      H      H      H   H   HHHHHHH   HHH     H        HHHH       HHHH        //
//      HH     HH      H      HH     H   H   H        H   H    H        H   HH    H    H       //
//      H H   H H     H H     H H    H   H   H       H     H   H        H     H    H           //
//      H  H H  H    H   H    H  H   H   H   H       H     H   H        H     H     HHH        //
//      H   H   H   HHHHHHH   H   H  H   H   HHHHH   H     H   H        H     H        H       //
//      H       H   H     H   H    H H   H   H       H     H   H        H     H         H      //
//      H       H  H       H  H     HH   H   H        H   H    H        H   HH    H   HH       //
//      H       H  H       H  H      H   H   H         HHH     HHHHHHH  HHHH       HHH         //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////


// Manifolds are objects which mimick the idea of a Riemannian manifold.
// A little bit of theory : a Riemannian manifold is a differentiable manifold
// endowed with a Riemannian metric. That is, for each point of the manifold
// an inner product is given on the tangent space to that manifold at that point.
// This induces a structure of metric space (a distance) on the tangent space and
// (through geodesics) on the manifold itself.


class ManiFEM::Manifold

{	public :

	// 'dim' is the dimension of the manifold
	// Might be immersed in R^n, with dim <= n. n is the size of the vectors received and returned
	// by 'project' and by 'rotate_vector_with_90_deg'.
	short int dim;

	MultiDimField * coord_field;
	FunctionOnMesh::Function * coord_func;

	// 'inner_prod_pointer' computes the inner product between two vectors belonging to
	// the space tangent to the manifold at a certain point. Used in 'inner_prod'.
	double (*inner_prod_pointer) ( void * data, const std::vector <double> & P, 
		const std::vector <double> & vec1, const std::vector <double> & vec2 );
	
	// 'distance_pointer' computes the distance bewteen two points of the manifold.
	// Note that this is not necessarily the geodesic distance.
	// It should probably be the Euclidean distance in the immersion space. Used in 'distance'.
	double (*distance_pointer) ( void * data, const std::vector <double> & P, const std::vector <double> & Q );

	// 'project_pointer' computes, for some point x in the immersion space, e.g. R^n,
	// a projected point belonging to the manifold. Used in 'project'.
	// It is assumed that x is not far away from the manifold, so a Newton method should do.
	std::vector <double> & (*project_pointer) ( void * data, const std::vector <double> & x );

	// 'rotate90_pointer' rotates 'vec' positively with 90º in the two-dimensional space
	// tangent to 'this' manifold at the given 'point'
	// (which means that the tangent space is oriented).
	// Plus, the resulting vector will have (Riemannian) norm equal to the argument 'norm'.
	// More precisely, the returned vector will be a norm-one vector times the argument 'norm'.
	// This is important since 'norm' may be negative, in which case we get a negative 90º rotation.
	// Used in 'rotate_vector_with_90_deg'.
	std::vector <double> & (*rotate_90_pointer) ( void * data, const std::vector <double> & point, 
		const std::vector <double> & vec, double desired_norm );

	// constructor
	
	inline Manifold ()
	{	Mesh::environment = this;  }
	
	inline Manifold ( short int d ) : Manifold ()
	{	this->dim = d;  }

	// 'inner_prod' computes the inner product between two vectors belonging to 
	// the space tangent to the manifold at a certain point.
	// The user should provide this function through 'inner_prod_pointer'.
	// This function may be invoked with PointerVectors instead of vector<double>,
	// based on the conversion operator defined in class PointerVector.
	// Not very efficient; is there another solution ?
	inline double inner_prod (const std::vector <double> & P,
	  const std::vector <double> & vec1, const std::vector <double> & vec2 )
	{	return inner_prod_pointer ( (void*) this, P, vec1, vec2 );  }

	// 'distance' computes the distance bewteen two points of the manifold.
	// Note that this is not necessarily the geodesic distance.
	// In some cases it could be the Euclidean distance in the immersion space,
	// but this depends on the user's intention.
	// The user should provide this function through 'distance_pointer'.
	// This function may be invoked with PointerVectors instead of vector<double>,
	// based on the conversion operator defined in class PointerVector.
	// Not very efficient; is there another solution ?

	inline double distance ( const std::vector <double> & P, const std::vector <double> & Q )
	{	return distance_pointer ( (void*) this, P, Q );  }

	// 'project' computes, for some point x in the immersion space, e.g. R^n,
	// a projected point belonging to the manifold.
	// The user should provide this function through 'project_pointer'.
	// It is assumed that x is not far away from the manifold, so a Newton method should do.
	// This function may be invoked with PointerVectors instead of vector<double>,
	// based on the conversion operator defined in class PointerVector.
	// Not very efficient; is there another solution ?

	inline std::vector <double> & project ( const std::vector <double> & x )
	{	return project_pointer ( (void*) this, x );  }

	// 'rotate_vector' rotates 'vec' positively with 90º in the two-dimensional space
	// tangent to 'this' manifold at the given 'point'
	// (which means that the tangent space is oriented).
	
	// Plus, the resulting vector will have (Riemannian) norm equal to the argument 'norm'.
	// More precisely, the returned vector will be a norm-one vector times the argument 'norm'.
	// This is important since 'norm' may be negative, in which case we get a negative 90º rotation.
	// The user should provide this function through 'rotate_90_pointer'.
	// This function may be invoked with PointerVectors instead of vector<double>,
	// based on the conversion operator defined in class PointerVector.
	// Not very efficient; is there another solution ?

	inline std::vector <double> & rotate_vector_with_90_deg
		( const std::vector <double> & point, const std::vector <double> & vec, double norm )
	{	assert ( dim == 2 );
		return rotate_90_pointer ( (void*) this, point, vec, norm );  };

	// factory functions

	static Manifold euclid ( short int d );
	Manifold & implicit ( FunctionOnMesh::Function & psi );
		
	inline void interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2 )
	{	MultiDimField & coord = * coord_field;
		coord.interpolate ( final, cell_1, frac_1, cell_2, frac_2 );
		PointerVector<double> & cf = coord(final);
		std::vector<double> cfp = project ( coord(final) );
		for ( size_t i = 0; i < cfp.size(); i++ ) cf[i] = cfp[i];       }
	
	inline void interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2,
	               Cell &cell_3, double frac_3                               )
	{	MultiDimField & coord = * coord_field;
		coord.interpolate ( final, cell_1, frac_1, cell_2, frac_2, cell_3, frac_3 );
		PointerVector<double> & cf = coord(final);
		std::vector<double> cfp = project ( coord(final) );
		for ( size_t i = 0; i < cfp.size(); i++ ) cf[i] = cfp[i];                    }
	
	inline void interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2,
	               Cell &cell_3, double frac_3, Cell &cell_4, double frac_4 )
	{	MultiDimField & coord = * coord_field;
		coord.interpolate ( final, cell_1, frac_1, cell_2, frac_2,
		                           cell_3, frac_3, cell_4, frac_4 );
		PointerVector<double> & cf = coord(final);
		std::vector<double> cfp = project ( coord(final) );
		for ( size_t i = 0; i < cfp.size(); i++ ) cf[i] = cfp[i];       }
	
	inline void interpolate
	( Cell &final, Cell &cell_1, double frac_1, Cell &cell_2, double frac_2,
	               Cell &cell_3, double frac_3, Cell &cell_4, double frac_4,
	               Cell &cell_5, double frac_5, Cell &cell_6, double frac_6 )
	{	MultiDimField & coord = * coord_field;
		coord.interpolate ( final, cell_1, frac_1, cell_2, frac_2,
		                           cell_3, frac_3, cell_4, frac_4,
		                           cell_5, frac_5, cell_6, frac_6 );
		PointerVector<double> & cf = coord(final);
		std::vector<double> cfp = project ( coord(final) );
		for ( size_t i = 0; i < cfp.size(); i++ ) cf[i] = cfp[i];       }
	
	// the following functions are used by the trivial manifold R^d (Euclidean space)
	static double trivial_inner_prod ( void*, const std::vector <double> & P,
		const std::vector <double> & vec1, const std::vector <double> & vec2 );
	static double trivial_distance ( void*, const std::vector <double> & P, const std::vector <double> & Q );
	static std::vector <double> & trivial_projection ( void*, const std::vector <double> & x );
	static std::vector <double> & trivial_rotate_90 ( void*, const std::vector <double> & point,
		const std::vector <double> & vec, double desired_norm );

	// the following functions are used by a submanifold defined through an implicit equation
	// (level set)
	static double implicit_mani_inner_prod ( void*, const std::vector <double> & P,
		const std::vector <double> & vec1, const std::vector <double> & vec2 );
	static double implicit_mani_distance
		( void*, const std::vector <double> & P, const std::vector <double> & Q );
	static std::vector <double> & implicit_mani_projection ( void*, const std::vector <double> & x );
	static std::vector <double> & implicit_mani_rotate_90 ( void*, const std::vector <double> & point,
		const std::vector <double> & vec, double desired_norm );
	
	// builds a FunctionOnMesh, with a MultiDimField associated
	// 's' may be "Lagrange degree one" for instance
	FunctionOnMesh::Function & coordinate_system ( std::string s );

}; // end of class ManiFEM::Manifold


class ManiFEM::ImplicitManifold : public Manifold

// a manifold defined by an implicit equation
// it's like a regular manifold plus a level function

{	public :
	
	Manifold * outer_space;

	// the projection will be done by means of the Newton method
	static const short int steps_for_Newton = 4;

	// '*level_function' should be an expression depending on the coordinates in '*outer_space'
	
	FunctionOnMesh::Function * level_function;
	FunctionOnMesh::Function deriv_lev_func;

	// constructor

	ImplicitManifold ( size_t n ) : Manifold(), deriv_lev_func(n)  { }
	
}; // end of class ManiFEM::ImplicitManifold



//////////////////////////////////////////////////////////////////////////////////
////////////////        variational problems     /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class ManiFEM::VariationalProblem

// an equality of the type " lhs == rhs "
// 'lhs' should be a bilinear form in unknown and test
// 'rhs' should be a linear form in test
// plus boundary conditions

{	public :

	// main equality
	Equality < FunctionOnMesh::combinIntegrals,
		FunctionOnMesh::combinIntegrals > equation;

	std::set < FunctionOnMesh::baseFunction* > space_coords, unknown, test;

	FiniteElement * finite_element;

	struct bdry_cond
	{	Equality < FunctionOnMesh::Function, FunctionOnMesh::Function > eq;
		Mesh * domain;                                                       };
	// boundary conditions
	std::list < bdry_cond * > conditions;

	static VariationalProblem * current_var_pb;
	static Mesh * current_mesh;

	// constructor
	inline VariationalProblem ( FunctionOnMesh::combinIntegrals * rhs,
															FunctionOnMesh::combinIntegrals * lhs )
		: equation ( rhs, lhs ) { }
	
	inline void prescribe_on ( Mesh & msh );
	// prepares the ground for declaring Dirichlet boundary conditions
	// followed by FunctionOnMesh::Function::operator==
	// defined after class ManiFEM::FiniteElement

	inline void set_finite_element ( FiniteElement & fe );
	// defined after class ManiFEM::FiniteElement

	void discretize ();

	void add_to_matrix ( FunctionOnMesh::Integral* );
	void add_to_vector ( FunctionOnMesh::Integral* );

}; // end of class ManiFEM::VariationalProblem
	
inline void FunctionOnMesh::Unknown::seek_content ( VariationalProblem * vp )
{	vp->unknown.insert (this);  }

inline void FunctionOnMesh::Test::seek_content (	VariationalProblem * vp )
{	vp->test.insert (this);  }

inline void FunctionOnMesh::Derivative::seek_content ( VariationalProblem * vp )
{	vp->space_coords.insert ( var );
	this->base->seek_content ( vp );    }

inline VariationalProblem	& FunctionOnMesh::combinIntegrals::operator==
	( FunctionOnMesh::combinIntegrals & rhs )
// factory function for variational problems
{	VariationalProblem * vp = new VariationalProblem ( this, &rhs );
	FunctionOnMesh::combinIntegrals * left = vp->equation.lhs;
	FunctionOnMesh::combinIntegrals * right = vp->equation.rhs;
	right->seek_content ( vp );
	assert ( vp->unknown.size() == 0 );  // no unknown in the right hand side
	assert ( vp->test.size() == 1 );  // test found in the right hand side
	left->seek_content ( vp );
	assert ( vp->unknown.size() == 1 );  // unknown found in the left hand side
	assert ( vp->test.size() == 1 );  // same test found in the left hand side
	return * vp;                                                                 }


//////////////////////////////////////////////////////////////////////////////////
//////////////////////      finite elements   ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace ManiFEM { namespace hidden {
		class FEwithMaster; class LagrangeFE;
} }

class ManiFEM::FiniteElement

// container for several types of finite elements

// This class is rather poor.
// Each finite element will have its own specific attributes.

{	public :

	// we use this class as a namespace :
	typedef hidden::FEwithMaster withMaster;
	typedef hidden::LagrangeFE Lagrange;

	static Cell * current_cell;

	std::list < FunctionOnMesh::baseFunction* > space_coords;
	std::list < FunctionOnMesh::baseFunction* > local_base;

	std::map < short int, Integrator * > integrator;

	// constructor

	inline FiniteElement ()
	{	// we make a copy of the list of space coordinates :
		std::vector<FunctionOnMesh::baseFunction*>::iterator it;
		for ( it = Mesh::environment->coord_func->components.begin();
					it != Mesh::environment->coord_func->components.end(); it++ )
			space_coords.push_back ( *it );
	}

	// factory function

	static FiniteElement::Lagrange & lagrange ( std::string );
	static FiniteElement::Lagrange & lagrange ( std::string, Mesh & msh );
	// the version without msh is called before building the mesh
	// reserves space for enumeration of future cells
	// the version with msh is called when the mesh already exists
	// re-sizes the size_t heap of existing cells for enumeration

	virtual void dock_on ( Cell & cll ) = 0;
	virtual void dock_on ( Cell & cll, std::string ) = 0;

	virtual inline size_t get_dof ( Cell & cll )
	{	assert (0);  }
	
	void set_integrator ( std::string, short int, std::string, short int, std::string );

	inline FunctionOnMesh::Function & new_coord ();
	
	inline virtual FunctionOnMesh::Integral * replace_space_coords
		( FunctionOnMesh::Integral * integr )
	{	return integr;  }
	
	inline virtual FunctionOnMesh::baseFunction * replace_space_coords
		( FunctionOnMesh::baseFunction * f )
	{	return f;  }
	// For elements which use no master, no need for replacing anything.
	// For elements using a master, see the method
	// FiniteElement::withMaster::replace_space_coords
	// declared below in class hidden::FEwithMaster and defined in Mesh.cpp

}; // end of class ManiFEM::FiniteElement

inline void VariationalProblem::set_finite_element ( FiniteElement & fe )
{	this->finite_element = & fe;  }

inline FunctionOnMesh::Function & FiniteElement::new_coord ()
{	FunctionOnMesh::Function *f = new FunctionOnMesh::Function (1);
	f->components[0] = new FunctionOnMesh::baseFunction;
	return *f;                                                        }

inline void FunctionOnMesh::prescribe_on ( Cell & cll )
// prepares the ground for setting values of a function
// followed by FunctionOnMesh::Function::operator==
{	VariationalProblem::current_var_pb = NULL;
	Mesh::current_cell_for_prescribe = &cll;         }

inline void VariationalProblem::prescribe_on ( Mesh & msh )
// declares Dirichlet boundary conditions
// followed by FunctionOnMesh::Function::operator==
{	FiniteElement::current_cell = NULL;
	VariationalProblem::current_var_pb = this;
	VariationalProblem::current_mesh = &msh;     }

inline void FunctionOnMesh::Function::operator== ( FunctionOnMesh::Function & f )
// used for imposing Dirichlet boundary conditions
{	assert ( VariationalProblem::current_var_pb != NULL );
	assert ( FiniteElement::current_cell == NULL );
	VariationalProblem::bdry_cond *bc = new VariationalProblem::bdry_cond;
	bc->eq.lhs = this;
	bc->eq.rhs = &f;
	bc->domain = VariationalProblem::current_mesh;
	VariationalProblem::current_var_pb->conditions.push_back(bc);          }

inline void FunctionOnMesh::Function::operator== ( double c )

// the equality operator between a function and a constant number
// can be used for two different purposes :
// either for imposing Dirichlet boundary conditions in a variational problem
// or for defining coordinates of a cell (e.g. a point)

{	if ( VariationalProblem::current_var_pb )
	{	(*this) == FunctionOnMesh::constant(c);  return;  }
	assert ( this->components.size() == 1 );
	*(this->components[0]) == c;                          }


inline void FunctionOnMesh::baseFunction::operator== ( double c )

// prescribe a value for a function on the current cell
	
// There are two possible situations.
// First : 'this' is a function with no field attached.
// So, we just prescribe a value for 'this' on the current cell.
// Second : there is a NumericField attached to 'this' function.
// In this case, we may want to set the value of that field.

{	assert ( Mesh::current_cell_for_prescribe != NULL );
	short int found_field = 0;
	std::list<RelationFunc*>::iterator it;
	RelationFuncField * copy_of_rel;
	for ( it = this->relations.begin(); it != this->relations.end(); it++ )
	{	if ( (*it)->type != "field" ) continue;
		RelationFuncField * rel = (RelationFuncField*) *it;
		AnyCell * ac = rel->cell;
		if ( ac->dim != Mesh::current_cell_for_prescribe->dim )
			continue;
		if ( ( ac->env_mesh == NULL ) && ( ac->env_cell == NULL ) )
			// universal cell
			{	found_field++; copy_of_rel = rel; continue;  }
		if ( ac->env_mesh )
			// we must check the current cell belongs to ac->env_mesh
		{	std::map < Mesh *, Cell::field_to_meshes > * cm =
				Mesh::current_cell_for_prescribe->meshes[ac->env_mesh->dim];
			if ( cm->find(ac->env_mesh) != cm->end() ) // yes, it belongs
			{	found_field++; copy_of_rel = rel; continue;  }                        }
	} // end for over this->relations
	assert ( ( found_field == 0 ) || ( found_field == 1 ) );
	if ( found_field )
	{	// simply set the value of that field on current cell
		OneDimField * fi = copy_of_rel->field;
		(*fi) ( * ( Mesh::current_cell_for_prescribe ) ) = c;
		return;                                      }
	// there is no field, just prescribe the value
	AnyCell * ac = new AnyCell ( Mesh::current_cell_for_prescribe );
	RelationFuncConst * rel = new RelationFuncConst;
	rel->val = c;
	rel->cell = ac;
	rel->func = this;
	this->relations.push_back (rel);                                           }
		
////////////////////////////////////////////////////////////////////////////////

class hidden::FEwithMaster : public FiniteElement

// also known as class FiniteElement::withMaster

{	public :

	std::list < FunctionOnMesh::baseFunction* > master_coords, composed_space_coords;

	FunctionOnMesh::Mapping * transf;

	// constructor

	inline FEwithMaster () : FiniteElement() { }

	virtual inline FunctionOnMesh::Integral * replace_space_coords
		( FunctionOnMesh::Integral * integr );

	virtual FunctionOnMesh::baseFunction * replace_space_coords
		( FunctionOnMesh::baseFunction * f );

}; // end of class hidden::FEwithMaster

inline FunctionOnMesh::Integral * FiniteElement::withMaster::replace_space_coords
( FunctionOnMesh::Integral * integr )
{	FunctionOnMesh::Integral * res = new FunctionOnMesh::Integral;
	res->domain = integr->domain;
	res->integrand = this->replace_space_coords ( integr->integrand );
	// calls replace_space_coords ( FunctionOnMesh::baseFunction * f )
	// (defined in Mesh.cpp)
	return res;                                                        }

//////////////////////////////////////////////////////////////////////////////////
////////////////        integrators        ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace ManiFEM { namespace hidden {
		class GaussIntegrator;
		double GaussQuadrature ( Integrator * integr, FunctionOnMesh::baseFunction & f );
} }

class ManiFEM::Integrator

// a container for Gauss quadrature (see class hidden::GaussIntegrator)
// and perhaps other integration methods (symbolic integration?)

{	public :

	// we use this class as a namespace
	typedef hidden::GaussIntegrator Gauss;
	
	FiniteElement * finite_element;
	
	std::list < FunctionOnMesh::baseFunction * > variables;

	double (*action) ( Integrator *, FunctionOnMesh::baseFunction & );
	// 'action' is strongly dependent on the type of integration we want
	// it could be Gauss quadrature with so many points
	// or it could be exact integration (through symbolic primitivation)

	// factory function

	static inline Integrator & gauss ( std::string type, std::string from = "main" );
	// sometimes this factory function is called directly from main,
	// sometimes it is called from FiniteElement::set_integrator
	// defined after class FiniteElement::Lagrange

	// operator()
	
	inline double operator() ( FunctionOnMesh::baseFunction & f )
	// When we receive just a function 'f' (as opposed to an integral)
	// it means that this call comes from a finite element which is already
	// docked on a cell, usually through a 'discretize' directive of some
	// variational formulation.. Thus, we just call 'action' to do the job.
	{	return (*action) ( this, f );  }
	// 'action' is strongly dependent on the type of integration we want
	
	inline double operator() ( FunctionOnMesh::Function & f )
	{	assert ( f.components.size() == 1 );
		return operator() ( * f.components[0] );  }
	
	//	inline double operator() ( FunctionOnMesh::combinIntegrals & in )
	// 'in' is a linear combination of integrals
	// (often just one integral wrapped as 'combinIntegrals')
	// we iterate over its terms and apply to each term the method
	// Integrator::operator(), which in turn calls 'action'
	//	{	std::list<pair<double,FunctionOnMesh::Integral*> >::iterator it;
	//	double res;
	//	for ( it = in.terms.begin(); it != in.terms.end(); it++ )
	//	res += it->first * this->operator() ( * ( it->second->integrand ) );
			// calls operator() ( FunctionOnMesh::baseFunction )
	//	return res;                                                 }

}; // end of class ManiFEM::Integrator


inline double FunctionOnMesh::baseFunction::integrate
( Mesh & msh, Integrator & integ )

// Here the integrator is called directly,
// not through a 'discretize' directive of some variational formulation.
// Thus, we need to sweep the entire mesh, to dock the finite element on each cell
// and all the stuff.

{	double res;
	FiniteElement & fe = * ( integ.finite_element );
	CellIterator it = msh.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
	for ( it.reset(); it.in_range(); it++ )
	// we add to 'res' the contribution of the cell '*it'
	{	Cell & cll = *it;
		fe.dock_on ( cll );
		// if a master element is used, we need to replace
		// each space_coord with the corresponding composed_space_coord
		// (if we use no master element, do nothing)
		FunctionOnMesh::baseFunction * f = fe.replace_space_coords ( this );
		// since this probably does not come from the 'discretize' directive of any
		// variational formulation, there should be no unknown function, no test function
		// also, we do not bother with delayed derivatives
		res += integ ( *f );                                                               }
	return res;                                                                             }


inline double FunctionOnMesh::Function::integrate
( Mesh & msh, Integrator & integ )
{	assert ( this->components.size() == 1 );
	return this->components[0]->integrate ( msh, integ ); } 

//////////////////////////////////////////////////////////////////////////////////

class hidden::GaussIntegrator : public Integrator

// also known as class Integrator::Gauss

// models a Gauss quadrature
// the class is flexible enough to model different quadratures
// on different space dimensions and with different number of points

// often, the integration is done on a master element by means of
// a coordinate transform, that's why we include the 'jacobian' attribute
// if we use no master element, the jacobian can be set to one or simply ingnored

{	public :

	FunctionOnMesh::baseFunction * jacobian;
  
	std::vector < double > gauss_weight;
	std::vector < std::vector < double > > gauss_point_coord;

};

//////////////////////////////////////////////////////////////////////////////////
////////////////   specific finite elements   ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


namespace ManiFEM { namespace hidden {
		class LagrangeFE_deg_1; class LagrangeFE_Q1;
		class LagrangeFE_deg_2; class LagrangeFE_Q2;
} }

class hidden::LagrangeFE : public FiniteElement::withMaster

// also known as class FiniteElement::Lagrange

{	public :

	// constructor

	inline LagrangeFE () : FiniteElement::withMaster() { }

	typedef hidden::LagrangeFE_deg_1 degree_1;
	typedef hidden::LagrangeFE_Q1 Q1;
	typedef hidden::LagrangeFE_deg_2 degree_2;
	typedef hidden::LagrangeFE_Q2 Q2;

}; // end of FiniteElement::Lagrange

//////////////////////////////////////////////////////////////////////////////////

class hidden::LagrangeFE_deg_1 : public FiniteElement::Lagrange

// also known as class FiniteElement::Lagrange::degree_1

// these are finite elements in which the base functions are polynomials of degree 1
// the degrees of freedom are the values of the function at vertices

{	public :

	// std::list < FunctionOnMesh::Function * > psi;

	// Degrees of freedom (functions in a base of the discretized Hilbert space)
	// must have an index (a global enumeration). For this particular type of finite element,
	// a convenient way to enumerate dofs is to enumerate vertices of the mesh.
	// 'pointer_into_heap' will be used to place that value in the heap of each vertex.
	// 'max_dof' will be used to keep track of the index of the next vertex.
	
	size_t pointer_into_heap, max_dof;

	std::list < Cell * > dofs;

	// constructor

	inline LagrangeFE_deg_1 () : FiniteElement::Lagrange() { }

	static inline void enumerate_new_vertex ( Cell * cll, void * data )
	// used to define the index of a newly built vertex
	{	assert ( cll->dim == 0 );
		assert ( cll->is_positive() );
		FiniteElement::Lagrange::degree_1 * fe =
			( FiniteElement::Lagrange::degree_1 * ) data;
		cll->size_t_heap[fe->pointer_into_heap] = fe->max_dof;
		fe->max_dof++;                                                    }

	virtual inline size_t get_dof ( Cell & cll )
	{	return cll.positive->size_t_heap.at(pointer_into_heap);  }   // debug
	
}; // end of FiniteElement::Lagrange::degree_1

//////////////////////////////////////////////////////////////////////////////////

class hidden::LagrangeFE_Q1 : public FiniteElement::Lagrange::degree_1

// also known as class FiniteElement::Lagrange::Q1

// Lagrange finite elements of degree 1 on quadrilaterals

{	public :

	// constructor

	inline LagrangeFE_Q1 () : FiniteElement::Lagrange::degree_1() { }

	virtual void dock_on ( Cell & cll );
	virtual void dock_on ( Cell & cll, std::string );
	
}; // end of FiniteElement::Lagrange::Q1
	
//////////////////////////////////////////////////////////////////////////////////

class hidden::LagrangeFE_deg_2 : public FiniteElement::Lagrange

// also known as class FiniteElement::Lagrange::degree_2

// these are finite elements in which the base functions are polynomials of degree 1
// the degrees of freedom are the values of the function at vertices

{	public :

	// std::list < FunctionOnMesh::Function * > psi;

	// degrees of freedom (functions in a base of the discretized Hilbert space)
	// must have an index (a global enumeration)
	// for this particular type of finite element, a convenient way to enumerate dofs is
	// to enumerate vertices and segments and quadrilaterals of the mesh
	// 'pointer_into_heap[0]' will be used to place that value in the heap of each vertex
	// 'pointer_into_heap[1]' will be used to place that value in the heap of each segment
	// 'pointer_into_heap[2]' will be used to place that value in the heap of each quadrilateral
	// 'max_dof' will be used to keep track of the index of the next cell
	
	size_t max_dof;

	std::vector < size_t > pointer_into_heap;

	std::list < Cell * > dofs;

	// constructor

	inline LagrangeFE_deg_2 () : FiniteElement::Lagrange()
	{ }

	inline LagrangeFE_deg_2 (size_t d) : FiniteElement::Lagrange(),
	                                     pointer_into_heap(d)
	// for triangles, d=2 : we only need to enumerate vertices and segments
	// for quadrilaterals, d=3 : we also need to enumerate quadrilaterals
	{ }

	static inline void enumerate_new_vertex ( Cell * cll, void * data )
	// used to define the index of a newly built vertex
	{	assert ( cll->dim == 0 );
		assert ( cll->is_positive() );
		FiniteElement::Lagrange::degree_2 * fe =
			( FiniteElement::Lagrange::degree_2 * ) data;
		cll->size_t_heap[fe->pointer_into_heap[0]] = fe->max_dof;
		fe->max_dof++;                                                    }

	static inline void enumerate_new_segment ( Cell * cll, void * data )
	// used to define the index of a newly built segment
	{	assert ( cll->dim == 1 );
		assert ( cll->is_positive() );
		FiniteElement::Lagrange::degree_2 * fe =
			( FiniteElement::Lagrange::degree_2 * ) data;
		cll->size_t_heap[fe->pointer_into_heap[1]] = fe->max_dof;
		fe->max_dof++;                                                    }

	virtual inline size_t get_dof ( Cell & cll )
	{	return cll.positive->size_t_heap.at(pointer_into_heap.at(cll.dim));  }  // debug
	
}; // end of FiniteElement::Lagrange::degree_2

//////////////////////////////////////////////////////////////////////////////////

class hidden::LagrangeFE_Q2 : public FiniteElement::Lagrange::degree_2

// also known as class FiniteElement::Lagrange::Q2

// Lagrange finite elements of degree 2 on quadrilaterals

{	public :

	// constructor

	inline LagrangeFE_Q2 () : FiniteElement::Lagrange::degree_2 (3)
	{ } // in the above, 3 means we enumerate cells of dimension 0, 1 and 2

	static inline void enumerate_new_quadr ( Cell * cll, void * data )
	// used to define the index of a newly built quadrilateral
	{	assert ( cll->dim == 2 );
		assert ( cll->is_positive() );
		FiniteElement::Lagrange::degree_2 * fe =
			( FiniteElement::Lagrange::degree_2 * ) data;
		cll->size_t_heap[fe->pointer_into_heap[2]] = fe->max_dof;
		fe->max_dof++;                                                    }

	virtual void dock_on ( Cell & cll );
	virtual void dock_on ( Cell & cll, std::string );
	
}; // end of FiniteElement::Lagrange::Q2
	

inline Integrator & Integrator::gauss ( std::string type, std::string from )

// if 'from' is missing, it is set by default to "main"
	
// sometimes this factory function is called directly from main,
// sometimes it is called from FiniteElement::set_integrator
// in the former case we need to create a finite element (if we want a master element)
// in the latter case the finite element already exists

{	assert ( type == "Q9");
  Integrator::Gauss * quadr = new Integrator::Gauss;
	std::vector < std::vector < double > > & coord = quadr->gauss_point_coord;
	std::vector < double > & weight = quadr->gauss_weight;
	size_t d = 2, n = 9;
	weight.resize(n);
	coord .resize(n);
	for ( size_t i = 0; i < n; i++ ) coord[i].resize(d);
	const double sqrt3over5 = sqrt(0.6);
	weight[0] = weight[8] = weight[2] = weight[6] = 25./81.;
	weight[1] = weight[3] = weight[5] = weight[7] = 40./81.;
	weight[4] = 64./81.;
	// {25./81., 40./81., 25./81., 40./81, 64./81., 40./81, 25./81., 40./81, 25./81.};
	coord[0][0] = coord[0][1] = coord[1][1] = coord[2][1] = coord[3][0] =
		coord[6][0] = -sqrt3over5;
	coord[2][0] = coord[5][0] = coord[6][1] = coord[7][1] = coord[8][0] =
		coord[8][1] = sqrt3over5;
	coord[1][0] = coord[3][1] = coord[4][0] = coord[4][1] = coord[5][1] =
		coord[7][0] = 0.;
	// {{-sqrt3over5,-sqrt3over5}, {0, -sqrt3over5}, {sqrt3over5, -sqrt3over5},
	//  {-sqrt3over5,0}, {0,0}, {sqrt3over5, 0}, {-sqrt3over5,sqrt3over5},
	//  {0,sqrt3over5}, {sqrt3over5, sqrt3over5}};
	quadr->action = & hidden::GaussQuadrature;
	if ( from == "main" )
	{	quadr->finite_element = & ( FiniteElement::lagrange ("Q1, no enum") );
		quadr->finite_element->integrator[2] = quadr;                           }
	return * quadr;                                                                     }
		



	
// trilogia qatsi koyaanisqatsi
