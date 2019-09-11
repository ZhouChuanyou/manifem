// manifem/global.cpp 2019.09.08

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "assert.h"
#include "Mesh.h"


// In Mesh::segment, the first two arguments are positive points.
// This is not consistent with Cell::segment, where the user
// must provide a negative point (the base) then a positive point (the tip).
// It is also inconsistent with other factory functions like Mesh::rectangle
// where we provide faces with orientation compatible with the orientation
// of the future mesh.
// However, we think it is easier for the user to build chains of segments like
//   auto & A = Cell::point(); auto & B = Cell::point();
//   auto & AB = Mesh::segment ( A, B, 10 );
// rather than
//   auto & A = Cell::point(); auto & B = Cell::point();
//   auto & AB = Mesh::segment ( A.reverse(), B, 10 );

Mesh & Mesh::segment ( Cell & a, Cell & b, size_t n ) // static

{	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );
	
	assert ( a.is_positive() );
	assert ( b.is_positive() );
	
	Mesh & malha = * ( new Mesh (1) );
	malha.hook["first point"] = (void*) &a;
	malha.hook["last point"] = (void*) &b;
	Cell * prev_point = &a;
	for ( size_t i=1; i < n; ++i )
	{	Cell & p = Cell::point();
		double frac = double(i)/double(n);
		space.interpolate ( p, a, 1.-frac, b, frac );
		Cell & seg = Cell::segment ( prev_point->reverse(), p );
		seg.add_to (malha); 
		prev_point = &p;                                         }
	Cell & seg = Cell::segment ( prev_point->reverse(), b );
	seg.add_to (malha);
	return malha;                                                         }


Mesh & Mesh::triangle ( Mesh & AB, Mesh & BC, Mesh & CA )  // static

{	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );

	// sides must be split in the same number of segments :
	size_t N = AB.number_of ( tag::cells, tag::of_max_dim );
	assert ( N == BC.number_of ( tag::cells, tag::of_max_dim ) );
	assert ( N == CA.number_of ( tag::cells, tag::of_max_dim ) );

	Cell * A_p, * B_p, * C_p;
	{ // just a block of code for hiding 'it'
	CellIterator it = AB.iter_over ( tag::vertices, tag::along );
	it.reset(); A_p = & (*it);
	} { // just a block of code for hiding 'it'
	CellIterator it = AB.iter_over ( tag::vertices, tag::along, tag::reverse );
	it.reset(); B_p = & (*it);
	} { // just a block of code for hiding 'it'
	CellIterator it = BC.iter_over ( tag::vertices, tag::along );
	it.reset(); assert ( B_p == & (*it) );
	} { // just a block of code for hiding 'it'
	CellIterator it = BC.iter_over ( tag::vertices, tag::along, tag::reverse );
	it.reset(); C_p = & (*it);
	} { // just a block of code for hiding 'it'
	CellIterator it = CA.iter_over ( tag::vertices, tag::along );
	it.reset(); assert ( C_p == & (*it) );
	} { // just a block of code for hiding 'it'
	CellIterator it = CA.iter_over ( tag::vertices, tag::along, tag::reverse );
	it.reset(); assert ( A_p == & (*it) );
	} // just a block of code for hiding 'it'

	// we keep a list of horizontal segments (parallel to AB)
	// useful for the next layer of triangles
	std::list <Cell*> ground, ceiling;
	{ // just a block of code for hiding 'it'
	CellIterator it = AB.iter_over ( tag::segments, tag::along );
	for ( it.reset(); it.in_range(); it++ ) ground.push_back ( & (*it) );
	} // just a block of code for hiding 'it'

	// we shall use six (pointers to) points, two on AB, two on BC, two on CA
	// like shadows of the point currently buing built
	Cell * P_AB_p, * Q_AB_p = A_p, * P_BC_p = B_p,
		* Q_BC_p, * P_CA_p, * Q_CA_p = A_p;

	// empty two-dimensional mesh :
	Mesh & tri_mesh = * ( new Mesh (2) );

	for ( short int i = 1; i <= N; i++ ) // "vertical" movement
	{	// advance one level upwards and slightly right (parallel to CA)
		{ // just a block of code for hiding 'seg'
		Cell & seg = AB.cell_in_front_of ( * Q_AB_p );
		Q_AB_p = & ( seg.tip() );
		} { // just a block of code for hiding 'seg'
		Cell & seg = BC.cell_in_front_of ( * P_BC_p );
		P_BC_p = & ( seg.tip() );
		} { // just a block of code for hiding 'seg'
		Cell & seg = CA.cell_behind ( * Q_CA_p );
		P_CA_p = Q_CA_p = & ( seg.base().reverse() );
		} // just a block of code for hiding 'seg'
		P_AB_p = A_p;  Q_BC_p = C_p;
		Cell * Q_AB_pc = Q_AB_p;  // keep a copy of Q_AB_p
		std::list<Cell*>::iterator it_ground = ground.begin();
		Cell * previous_seg_p, * ground_ver_p;
		// build the first triangle on this layer
		{ // just a block of code for hiding variables
		Cell & AS = * ( * it_ground );
		Cell & S = AS.tip();
		ground_ver_p = &S;
		Cell & ST = Cell::segment ( S.reverse(), * P_CA_p );
		previous_seg_p = & ST;
		Cell & tri = Cell::triangle ( AS, ST, CA.cell_in_front_of(*P_CA_p) );
		tri.add_to(tri_mesh);
		} // just a block of code for hiding variables
		Cell * previous_ver_p = Q_CA_p;
		ceiling.clear();
		for ( short int j = i+1; j <= N; j++ ) // "horizontal" movement
		{	// advance one step horizontally (parallel to AB)
			{ // just a block of code for hiding 'seg'
			Cell & seg = AB.cell_in_front_of ( * P_AB_p );
			P_AB_p = & ( seg.tip() );
			} { // just a block of code for hiding 'seg'
			Cell & seg = AB.cell_in_front_of ( * Q_AB_pc );
			Q_AB_pc = & ( seg.tip() );
			} { // just a block of code for hiding 'seg'
			Cell & seg = BC.cell_behind ( * Q_BC_p );
			Q_BC_p = & ( seg.base().reverse() );
			} { // just a block of code for hiding 'seg'
			Cell & seg = CA.cell_behind ( * P_CA_p );
			P_CA_p = & ( seg.base().reverse() );
			} // just a block of code for hiding 'seg'
			Cell * new_ver_p;
			if ( j == N ) new_ver_p = P_BC_p;
			else
			{	// we prepare for building a new point and we need fractions
				// distance to AB : i
				// distance to BC : N-j
				// distance to CA : j-i
				double frac_AB = 1. / double(i),
				       frac_BC = 1. / double(N-j),
				       frac_CA = 1. / double(j-i);
				double s = frac_AB + frac_BC + frac_CA;  s *= 2.;
				frac_AB /= s;  frac_BC /= s;  frac_CA /= s;
				new_ver_p = & ( Cell::point() );
				space.interpolate ( *new_ver_p, *P_AB_p, frac_AB, *Q_AB_pc, frac_AB,
														            *P_BC_p, frac_BC, *Q_BC_p,  frac_BC,
														            *P_CA_p, frac_CA, *Q_CA_p,  frac_CA  );  }
			Cell & new_seg = Cell::segment ( ground_ver_p->reverse(), *new_ver_p );
			Cell & horizontal_seg = Cell::segment ( new_ver_p->reverse(), *previous_ver_p );
			Cell & tri_1 = Cell::triangle ( previous_seg_p->reverse(), new_seg, horizontal_seg );
			tri_1.add_to ( tri_mesh );
			it_ground++; assert ( it_ground != ground.end() );
			Cell & ground_seg = * (*it_ground);
			ground_ver_p = & ( ground_seg.tip() );
			previous_seg_p = & ( Cell::segment ( ground_ver_p->reverse(), *new_ver_p ) );
			Cell & tri_2 = Cell::triangle ( ground_seg, *previous_seg_p, new_seg.reverse() );
			tri_2.add_to ( tri_mesh );
			previous_ver_p = new_ver_p;
			// add horizontal_seg.reverse() to future ground
			ceiling.push_back ( & ( horizontal_seg.reverse() ) );                              	  }
		ground = ceiling;                                                                         }

	return tri_mesh;
	
} // end of Mesh::triangle

	
Mesh & Mesh::rectangle  // static
( Mesh & south, Mesh & east, Mesh & north, Mesh & west, const tag::WithTriangles & wt )

// 'wt' defaults to 'tag::not_with_triangles',
// which means 'cut_rectangles_in_half' defaults to 'false'

{	bool cut_rectangles_in_half = ( wt == tag::with_triangles );

	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );

	// recover the corners from the sides
	Cell *SE, *NE, *SW, *NW;
	size_t N_horiz, N_vert;
	// we get the first and last points of each side using freshly reset iterators
	// as explained in the manual, at the end of paragraph 2.5
	{ // just to keep it1 and it2 local
	CellIterator it1 = south.iter_over ( tag::vertices, tag::along );
	it1.reset(); SW = &(*it1);
	CellIterator it2 = south.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); SE = &(*it2);
	N_horiz = south.number_of ( tag::cells, tag::of_max_dim );
	} { // just to keep it1 and it2 local
	CellIterator it1 = north.iter_over ( tag::vertices, tag::along);
	it1.reset(); NE = &(*it1);
	CellIterator it2 = north.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); NW = &(*it2);
	assert ( N_horiz == north.number_of ( tag::cells, tag::of_max_dim ) );
	} { // just to keep it1 and it2 local
	CellIterator it1 = east.iter_over ( tag::vertices, tag::along );
	it1.reset(); assert ( SE == &(*it1) );
	CellIterator it2 = east.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); assert ( NE == &(*it2) );
	N_vert = east.number_of ( tag::cells, tag::of_max_dim );
	} { // just to keep it1 and it2 local
	CellIterator it1 = west.iter_over ( tag::vertices, tag::along );
	it1.reset(); assert ( NW == &(*it1) );
	CellIterator it2 = west.iter_over ( tag::vertices, tag::along, tag::reverse);
	it2.reset(); assert ( SW == &(*it2) );
	assert ( N_vert == west.number_of ( tag::cells, tag::of_max_dim ) );
	} // just to keep it1 and it2 local
		
	// prepare the mesh
	// empty two-dimensional mesh :
	Mesh & rectangle = * ( new Mesh (2) );
	std::map <std::string, Mesh*> * sides_p = new std::map <std::string, Mesh*>;
	std::map <std::string, Mesh*> & sides = * sides_p;
	sides["south"] = &south; sides["east"] = &east;
	sides["north"] = &north; sides["west"] = &west;
	rectangle.hook["sides"] = (void*) ( sides_p );
	std::map <std::string, Cell*> * corners_p = new std::map <std::string, Cell*>;
	std::map <std::string, Cell*> & corners = * corners_p;
	corners["SW"] = SW; corners["SE"] = SE; corners["NE"] = NE; corners["NW"] = NW;
	rectangle.hook["corners"] = (void*) ( corners_p );
	
	// prepare horizon
	std::list <Cell*> horizon;
	{ // just a block of code for hiding 'it'
	CellIterator it = south.iter_over ( tag::segments, tag::along );
	for ( it.reset(); it.in_range(); it++ )
	{	Cell & seg = *it;
		horizon.push_back ( &seg );  }
	} // just a block of code for hiding 'it'

	// start mesh generation
	CellIterator it_east = east.iter_over ( tag::vertices, tag::along );
	CellIterator it_west = west.iter_over ( tag::vertices, tag::along, tag::reverse );
	CellIterator it_south = south.iter_over ( tag::vertices, tag::along );
	CellIterator it_north = north.iter_over ( tag::vertices, tag::along, tag::reverse );
	it_east.reset(); it_east++;
  it_west.reset(); it_west++;
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list <Cell*>::iterator it = horizon.begin();
		Cell & A = (*it)->base().reverse();
		Cell * s4_p = west.cell_behind ( A, tag::may_not_exist );
		assert ( s4_p != NULL );
		Cell * D_p = & (  s4_p->base().reverse() );
		Cell & ver_east = *it_east;
		Cell & ver_west = *it_west;
		double frac_N = double(i) / double(N_vert),
					alpha = frac_N * (1-frac_N);
		alpha = alpha*alpha*alpha;
		it_south.reset(); it_south++;
		it_north.reset(); it_north++;
		for ( size_t j = 1; j < N_horiz; j++ )
		{	Cell * s1_p = *it;  // 'it' points into the 'horizon' list
			Cell & B = s1_p->tip();
			Cell & C = Cell::point();
			Cell & ver_south = *it_south;
			Cell & ver_north = *it_north;
			double frac_E = double(j) / double(N_horiz),
				beta = frac_E * (1-frac_E);
			beta = beta*beta*beta;
			double sum = alpha + beta,
				aa = alpha/sum,  bb = beta/sum;
			space.interpolate ( C, ver_south, bb*(1-frac_N), ver_east, aa*frac_E,     
			                    ver_north, bb*frac_N,     ver_west, aa*(1-frac_E));
			Cell & s2 = Cell::segment ( B.reverse(), C );
			Cell & s3 = Cell::segment ( C.reverse(), *D_p );
			if ( cut_rectangles_in_half )
			{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
				Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
				Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
				T1.add_to (rectangle);
				T2.add_to (rectangle);                                 }
			else // with quadrilaterals
			{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p );
				Q.add_to (rectangle);                                 }
			// 's3' is on the ceiling, we keep it in the 'horizon' list
			// it will be on the ground when we build the next layer of cells
			*it = & ( s3.reverse() ); // 'it' points into the 'horizon' list
			it++;
			D_p = &C;
			s4_p = & ( s2.reverse() );
			it_south++;  it_north++;
		} // end of for j
		Cell * s1_p = *it;
		Cell & B = s1_p->tip();
		Cell & s2 = east.cell_in_front_of (B);
		Cell & C = s2.tip();
		Cell & s3 = Cell::segment ( C.reverse(), *D_p );
		if ( cut_rectangles_in_half )
		{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
			Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
			Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
			T1.add_to (rectangle);
			T2.add_to (rectangle);                                 }
		else // with quadrilaterals
		{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p );
			Q.add_to (rectangle);                                 }
		*it = & ( s3.reverse() );
		it_east++;  it_west++;
	} // end of for i
	std::list <Cell*>::iterator it = horizon.begin();
	Cell * s4_p = west.cell_in_front_of ( *NW, tag::may_not_exist );
	assert ( s4_p != NULL );
	Cell * D_p = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	//Cell & s1 = south.cell_in_front_of (A);
		Cell * s1_p = *it;
		Cell & B = s1_p->tip();
		Cell & s3 = north.cell_behind ( *D_p );
		Cell & C = s3.base().reverse();
		Cell & s2 = Cell::segment ( B.reverse(), C );
		if ( cut_rectangles_in_half )
		{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
			Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
			Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
			T1.add_to (rectangle);
			T2.add_to (rectangle);                                 }
		else // with quadrilaterals
		{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p );
			Q.add_to (rectangle);                                 }
		it++;
		D_p = &C;
		s4_p = & ( s2.reverse () );                             }
	Cell * s1_p = *it;
	Cell & B = s1_p->tip();
	Cell & s2 = east.cell_in_front_of (B);
	Cell & C = s2.tip();
	assert ( &C == NE );
	Cell & s3 = north.cell_behind ( *D_p );
	if ( cut_rectangles_in_half )
	{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
		Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
		Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
		T1.add_to (rectangle);
		T2.add_to (rectangle);                                 }
	else // with quadrilaterals
	{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p);
		Q.add_to (rectangle);                                 }
	return rectangle;

} // end of Mesh::rectangle


Mesh & Mesh::join ( const std::list<Mesh*> & l ) // static

{	// compute the dimension
	std::list<Mesh*>::const_iterator it = l.begin();
	assert ( it != l.end() );
	Mesh & first = *(*it);
	short int dim = first.dim;
	bool compute_first_last = ( dim == 1 );
	for ( it++; it != l.end(); it++ )
	{	Mesh & other = *(*it);
		assert ( dim == other.dim );  }
	// prepare the mesh
	Mesh & result = * ( new Mesh (dim) );
	// sweep the list of meshes
	for ( it = l.begin(); it != l.end(); it++ )
	{	Mesh & m = *(*it);
		// sweep all cells of m
		CellIterator itt = m.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( itt.reset(); itt.in_range(); itt++ )
		{	Cell & cll = *itt;
			cll.add_to ( result );  }
		compute_first_last = compute_first_last
			&& ( m.hook.find("first point") != m.hook.end() )
			&& ( m.hook.find("last point") != m.hook.end() );                        	}
	if ( compute_first_last )
	{	std::set < void* > set_of_first, set_of_last;
		for ( it = l.begin(); it != l.end(); it++ )
		{	Mesh & m = *(*it);
			set_of_first.insert ( m.hook["first point"] );
			set_of_last.insert ( m.hook["last point"] );    }
		std::set<void*>::iterator it_first = set_of_first.begin(),
		                     it_last  = set_of_last.begin();
		short int found_first = 0, found_last = 0;
		for ( ; it_first != set_of_first.end(); it_first++, it_last++ )
		{	assert ( it_last != set_of_last.end() );
			if ( set_of_first.find(*it_last) == set_of_first.end() )
			{	found_last++; result.hook["last point"] = *it_last;  }
			if ( set_of_last.find(*it_first) == set_of_last.end() )
			{	found_first++; result.hook["first point"] = *it_first;  }		}
		assert ( it_last == set_of_last.end() );                           }
	return result;                                                                   }


Cell & Cell::cartesian_product_orient ( Cell & cell1, Cell & cell2,  // static
	std::map < Cell*, std::map < Cell*, std::pair < Cell*, bool > > > & cartesian, bool revert )

// receives two oriented cells, finds their positive counterparts,
// builds their cartesian product and then gives it the correct orientation
	
{	Cell * cell1_pos = cell1.positive;
	Cell * cell2_pos = cell2.positive;

	if ( cell1_pos != (&cell1) ) revert = not revert;
	if ( cell2_pos != (&cell2) ) revert = not revert;

	// recover_cartesian_product ( cell1_pos, cell2_pos, cartesian, revert );
	{ // just a block of code for hiding it1, p
	Cell * p = NULL;
	std::map<Cell*,std::map<Cell*,std::pair<Cell*,bool>>>::iterator it1 = cartesian.find (cell1_pos);
	if ( it1 != cartesian.end() )
	{	std::map<Cell*,std::pair<Cell*,bool>> & cc1 = it1->second; // cartesian[cell1_pos]
		std::map<Cell*,std::pair<Cell*,bool>>::iterator it2 = cc1.find(cell2_pos);
		if ( it2 != cc1.end() )
		{	if ( it2->second.second ) revert = not revert;
			// we are changing the value of revert, but it's OK to do that
			// because we will surely return at the end of this block
			if ( revert ) p = & ( it2->second.first->reverse() );
			else p = it2->second.first;                           }   }
	if ( p != NULL ) return *p;
	} // just a block of code for hiding it1, p

	// we give a special treatment to points since they have no boundary, and have coordinates
	
	if ( ( cell1.dim == 0 ) and ( cell2.dim == 0 ) )
	{	Cell & new_point = Cell::point();
		assert ( not revert );
		assert ( cell1.is_positive() );
		assert ( cell2.is_positive() );
		FunctionOnMesh::Function & xyz = * ( Mesh::environment->coord_func );
		FunctionOnMesh::Function & x = xyz[0], & y = xyz[1], & z = xyz[2];
		x == x(cell1); y == y(cell1); z == z(cell2);
		cartesian[cell1_pos][cell2_pos] = std::pair ( &new_point, revert );
		return new_point;                                                       }

	// segments must also be given a special treatment because they have
	// special attributes (base and tip) which do not exist for other cells

	if ( cell1.dim + cell2.dim == 1 )
	{
	Cell * first, * second;
	if ( ( cell1.dim == 1 ) and ( cell2.dim == 0 ) )
	{	first = & ( Cell::cartesian_product_orient
		    ( cell1_pos->base().reverse(), *cell2_pos, cartesian, false ) );
		second = & ( Cell::cartesian_product_orient
		    ( cell1_pos->tip(), *cell2_pos, cartesian, false ) );
		// 'false' means 'do not revert'
		assert ( first->is_positive() );
		assert ( second->is_positive() );                                     }
	if ( ( cell1.dim == 0 ) and ( cell2.dim == 1 ) )
	{	first = & ( Cell::cartesian_product_orient
		      ( *cell1_pos, cell2_pos->base().reverse(), cartesian, false ) );
		second = & ( Cell::cartesian_product_orient
		      ( *cell1_pos, cell2_pos->tip(), cartesian, false ) );
		// 'false' means 'do not revert'
		assert ( first->is_positive() );
		assert ( second->is_positive() );                                        }
	Cell * prod_cell;
	if ( revert )
		prod_cell = & ( Cell::segment ( second->reverse(), *first ) );
	else
		prod_cell = & ( Cell::segment ( first->reverse(), *second ) );
	cartesian[cell1_pos][cell2_pos] = std::pair ( prod_cell, revert );
	return *prod_cell;
	} // end if

	// dimension higher than 1

	Cell & prod_cell = * ( new Cell ( cell1.dim + cell2.dim ) );
	// we build the boundary, initially empty :
	Mesh::build_boundary ( prod_cell );

	// now we add faces to the boundary, one by one
	// "vertical" faces : faces of cell1 multiplied by cell2 itself
	if ( cell1.dim > 0 )
	{	CellIterator it1 = cell1_pos->boundary().iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( it1.reset(); it1.in_range(); it1++ )
		{	Cell & face = Cell::cartesian_product_orient ( *it1, *cell2_pos, cartesian, revert );
			assert ( face.dim == prod_cell.dim - 1 );
			face.glue_on_bdry_of ( prod_cell );                                                    }       }
	// "horizontal" faces : cell1 multiplied by faces of cell2
	if ( cell2.dim > 0 )
	{	bool local_revert = revert;
		if ( cell1.dim %2 ) local_revert = not local_revert;
		// the above works (for the cases I've tested), although I'm not sure why;
		// I could use some help from experts in algebraic topology or similar domains
		CellIterator it2 = cell2_pos->boundary().iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		for ( it2.reset(); it2.in_range(); it2++ )
		{	Cell & face = Cell::cartesian_product_orient ( *cell1_pos, *it2, cartesian, local_revert );
			assert ( face.dim == prod_cell.dim - 1 );
			face.glue_on_bdry_of ( prod_cell );                                                          }  }

	cartesian[cell1_pos][cell2_pos] = std::pair ( &prod_cell, revert );
	return prod_cell;

} // end of Cell::cartesian_product_orient


Mesh & Mesh::cartesian_product ( Mesh & mesh1, Mesh & mesh2 ) // static

{	// revert if mesh2.dim is odd
	bool revert = ( mesh2.dim % 2 );

	// the above works (for the cases I've tested), although I'm not sure why;
	// I could use some help from experts in algebraic topology or similar domains
	
	Mesh * product_mesh_p = new Mesh ( mesh1.dim + mesh2.dim );	 
	Mesh & product_mesh = *product_mesh_p;

	std::map < Cell*, std::map < Cell*, std::pair < Cell*, bool > > >  cartesian;

	CellIterator it1 = mesh1.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );	
	CellIterator it2 = mesh2.iter_over ( tag::cells, tag::of_max_dim, tag::oriented );

	for(it1.reset(); it1.in_range(); it1++)
	for(it2.reset(); it2.in_range(); it2++)
	{	Cell & cell1 = *it1;
		Cell & cell2 = *it2;
		Cell & prod_cell = Cell::cartesian_product_orient ( cell1, cell2, cartesian, false );
		// 'false' means 'do not revert'
		assert ( prod_cell.dim == product_mesh.dim );
		prod_cell.add_to ( product_mesh );                                                     }

	return product_mesh;

} // end of Mesh::cartesian_product


void Mesh::draw_ps ( std::string f )

// cabecalho para o seguinte comando funcionar :
// gs -q -dNOPAUSE -dBATCH -sDEVICE=png16 -sOUTPUTFILE=output.png input.eps
	
// %!PS-Adobe-3.0 EPSF-3.0
// %%Title:                     malha
// %%BoundingBox:  0 0  520   270
// %%EndComments
// %%BeginSetup
// << /PageSize [520 270] >> setpagedevice
// %%EndSetup
	
{	assert ( Mesh::environment != NULL );
	NumericField &coord = * ( Mesh::environment->coord_field );

	assert (coord.size() == 2);
	// we split 'coord' into its components
	OneDimField &x = coord[0], &y = coord[1];
	
	double xmin, xmax, ymin, ymax, maxside;
	
	{ // just a block for hiding variables
	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
	it.reset(); assert(it.in_range());
	Cell & Vfirst = *it;
	xmin = xmax = x(Vfirst);
	ymin = ymax = y(Vfirst);	
	for ( it++ ; it.in_range(); it++ )
	{ Cell & V = *it; 
		double xV = x(V), yV = y(V);
		if ( xV < xmin ) xmin = xV;
	  if ( xV > xmax ) xmax = xV;
	  if ( yV < ymin ) ymin = yV;
	  if ( yV > ymax ) ymax = yV;      }
	} // just a block for hiding variables 
	if ( xmax-xmin < ymax-ymin ) maxside = ymax-ymin;
	else maxside = xmax-xmin;
	double border = 0.02*maxside;
	double scale_factor = 500/maxside;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;

	std::ofstream file_ps (f);
	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     malha" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(ymax-ymin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(ymax-ymin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;
	
	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;

	file_ps << "gsave " << 1.5 / scale_factor << " setlinewidth" << std::endl;
	
	{ // just a block for hiding variables
	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 1, tag::not_oriented );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell & seg = *it;
		Cell & base = seg.base().reverse();
		Cell & tip  = seg.tip();
		file_ps << x(base) << " " << y(base) << " moveto" << std::endl;
		file_ps << x(tip) << " " << y(tip) << " lineto stroke" << std::endl;  }
	} // just a block for hiding variables
																					
#ifndef NDEBUG
	file_ps << "/Courier findfont 0.2 scalefont setfont" << std::endl;
	{ // just a block for hiding variables
	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell & p = *it;
		if ( p.name() == "" ) continue;
		file_ps << x(p) << " " << y(p) << " moveto (" << p.name() << ") show" << std::endl;  }
	} // just a block for hiding variables
#endif

	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	if ( ! file_ps.good() )
	{	std::cerr << "error writing postscript file" << std::endl;
		exit (1);                                    }

} // end of Mesh::draw_ps


void Mesh::export_msh ( std::string f, std::map<Cell*,size_t> & ver_numbering )

{	assert ( Mesh::environment != NULL );
	NumericField &coord = * ( Mesh::environment->coord_field );

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;
	
  file_msh << "$Nodes" << std::endl << this->number_of(tag::cells,tag::of_dim,0) << std::endl;

	{ // just to make variables local : it, counter, x, y
	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
	size_t counter = 0;
	OneDimField & x = coord[0], & y = coord[1];
	if (coord.size() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	++counter;
			Cell & p = *it;
			ver_numbering [&p] = counter;
			file_msh << counter << " " << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	if ( coord.size() == 3 )
	{	OneDimField &z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	++counter;
			Cell & p = *it;
			ver_numbering [&p] = counter;
			file_msh << counter << " " << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just to make variables local : it, counter, x, y

	file_msh << "$Elements" << std::endl;
	file_msh << this->number_of ( tag::cells, tag::of_max_dim ) << std::endl;

	if ( this->dim == 2 )
	{	CellIterator it = this->iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell & elem = *it;
			if ( elem.boundary().number_of ( tag::cells, tag::of_max_dim ) == 3 ) // a triangle
				file_msh << counter << " 2 0 ";
			else // a quadrilateral
			{	assert ( elem.boundary().number_of ( tag::cells, tag::of_max_dim ) == 4 );
				file_msh << counter << " 3 0 ";                            }
			CellIterator itt = elem.boundary().iter_over ( tag::vertices, tag::around );
			for ( itt.reset(); itt.in_range(); ++itt )
			{	Cell & p = *itt;
				file_msh << ver_numbering [&p] << " ";   }
			file_msh << std::endl;                                                            }  }
	else
	{	assert ( this->dim == 3);
		CellIterator it = this->iter_over ( tag::cells, tag::of_max_dim, tag::oriented );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell & elem = *it;
			size_t n_faces = elem.boundary().number_of ( tag::cells, tag::of_max_dim );
			if ( n_faces == 4 ) // a tetrahedron
				file_msh << counter << " 4 0 ";  // to finish !
			else if ( n_faces == 6 )
			{	// 3d parallelogram = 8-node hexaedron = cube
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 5 0 ";
				CellIterator itt = elem.boundary().iter_over
					( tag::cells, tag::of_max_dim, tag::oriented );
				itt.reset();
				Cell & back = *itt; // square face behind the cube
				// back is 0321 in gmsh's documentation
				assert ( back.boundary().number_of ( tag::cells, tag::of_max_dim ) == 4 );
				CellIterator itv = back.boundary().iter_over
					( tag::vertices, tag::around, tag::reverse );
				// reverse because we want the vertices ordered as 0, 1, 2, 3
				itv.reset();  Cell & ver_0 = *itv;
				Cell & seg_03 = back.boundary().cell_in_front_of(ver_0);
				for ( ; itv.in_range(); ++itv )
				{	Cell & p = *itv;
					file_msh << ver_numbering [&p] << " ";   }
				Cell & left_wall = elem.boundary().cell_in_front_of(seg_03); // square face on the left
				// left_wall is 0473 in gmsh's documentation
				assert ( left_wall.boundary().number_of ( tag::cells, tag::of_max_dim ) == 4 );
				Cell & seg_04 = left_wall.boundary().cell_in_front_of(ver_0);
				Cell & ver_4 = seg_04.tip();
				Cell & seg_47 = left_wall.boundary().cell_in_front_of(ver_4);
				Cell & front = elem.boundary().cell_in_front_of(seg_47); // square face in front
				// front is 4567 in gmsh's documentation
				assert ( front.boundary().number_of ( tag::cells, tag::of_max_dim ) == 4 );
				CellIterator itvv = front.boundary().iter_over ( tag::vertices, tag::around );
				itvv.reset(ver_4);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell & p = *itvv;
					file_msh << ver_numbering [&p] << " ";   }                                         }
			else
			{	assert( n_faces == 5 );
				// triangular prism = 6-node prism
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 6 0 ";
				CellIterator itt = elem.boundary().iter_over
					( tag::cells, tag::of_max_dim, tag::oriented );
				short int n_tri = 0, n_rect = 0;
				Cell* base_p;
				for( itt.reset(); itt.in_range(); ++itt )
				{	Cell & face = *itt; // every face elem
			    size_t n_edges = face.boundary().number_of ( tag::cells, tag::of_max_dim );
					if ( n_edges == 3 )  { n_tri++;   base_p = &face;           }
					else                 { n_rect++;  assert ( n_edges == 4 );  }   }
				assert ( n_tri == 2 );  assert ( n_rect == 3 );
				// base is 021 in gmsh's documentation
				assert ( base_p->boundary().number_of ( tag::cells, tag::of_max_dim ) == 3 );
				CellIterator itv = base_p->boundary().iter_over
					( tag::vertices, tag::around, tag::reverse );
				// we set it reverse because we want the vertices ordered as 0, 1, 2
				itv.reset();
				Cell & ver_0 = *itv;
				Cell & seg_02 = base_p->boundary().cell_in_front_of(ver_0);
				for (  ; itv.in_range(); ++itv )
				{	Cell & p = *itv;
					file_msh << ver_numbering [&p] << " ";   }
				Cell & right_wall = elem.boundary().cell_in_front_of(seg_02);
				// right_wall is 0352 in gmsh's documentation
				assert ( right_wall.boundary().number_of ( tag::cells, tag::of_max_dim ) == 4 );
				Cell & seg_03 = right_wall.boundary().cell_in_front_of(ver_0);
				Cell & ver_3 = seg_03.tip();
				Cell & seg_35 = right_wall.boundary().cell_in_front_of(ver_3);
				Cell & roof = elem.boundary().cell_in_front_of(seg_35); //
				// roof is 345 in gmsh's documentation
				assert ( roof.boundary().number_of ( tag::cells, tag::of_max_dim ) == 3 );
				CellIterator itvv = roof.boundary().iter_over ( tag::vertices, tag::around );
				itvv.reset(ver_3);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell & p = *itvv;
					file_msh << ver_numbering [&p] << " ";   }                                      }
			file_msh << std::endl;                                                                } }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f )
	
// the numbering of vertices is produced on-the-fly

{	std::map < Cell*, size_t > numbering;

	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++ )
	{	++counter;  Cell & p = *it;  numbering [&p] = counter;  }

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f, FiniteElement & fe )

// the numbering of vertices is given by the degrees of freedom provided by 'fe'

{
	std::map < Cell*, size_t > numbering;

	CellIterator it = this->iter_over ( tag::cells, tag::of_dim, 0, tag::not_oriented );
	for ( it.reset() ; it.in_range(); it++ )
	{	Cell & p = *it;  numbering [&p] = fe.get_dof(p)+1;  }

	this->export_msh ( f, numbering );

} // end of Mesh::export_msh


Cell & Cell::split ( const tag::InTwoSegments &, const tag::CellIsSegment & )

// split 'this' segment in two smaller segments
// after split, 'this' will be the first half of the segment
// return the second half

{	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );
	
	Cell & A = this->base().reverse();
	Cell & B = this->tip();
	
	Cell & middle = Cell::point();
	space.interpolate ( middle, A, 0.5, B, 0.5 );
	B. cut_from_bdry_of ( *this );
	middle.glue_on_bdry_of ( *this );
	Cell & AB2 = Cell::segment ( middle.reverse(), B );
	
	// we must introduce AB2 to all meshes above 'this'
	AB2.is_part_of ( *this );

	return AB2;
}

		
Cell & Cell::split ( const tag::InTwoTriangles &, const tag::At &,
                   Cell & ver, const tag::CellIsRectangle &     )

// split 'this' rectangle along a diagonal
// which diagonal ? the one passing through vertex 'ver'
// returns the newly created segment, beginning at 'ver'

// split_from_bdry and glue_to_bdry deal with the reverse on their own
// add_to does too, I think (to check)

{	Mesh & bdry = this->boundary();
	// let's make sure 'this' is a rectangle (quadrilateral, actually)
	assert ( this->dim == 2 );
	assert ( bdry.number_of(tag::cells,tag::of_max_dim) == 4 );

	// unglue two segments from the boundary
	Cell & seg1 = bdry.cell_in_front_of (ver);
	seg1.cut_from_bdry_of ( * this );
	// 'this' remains incomplete for now, open to the exterior
	Cell & ver1 = seg1.tip();
	Cell & seg2 = bdry.cell_in_front_of (ver1);
	seg2.cut_from_bdry_of ( * this );
	// 'this' becomes widely open, only two walls left
	Cell & ver2 = seg2.tip();

	// now close 'this' back
	Cell & new_seg = Cell::segment ( ver.reverse(), ver2 );
	new_seg.glue_on_bdry_of ( * this );
	// 'this' is closed now, has three walls (has become a triangle)

	// now build a new triangle, to take account of the missing half
	// of the original 'this'
	Cell & new_tri = Cell::triangle ( seg1, seg2, new_seg.reverse() );

	// if 'this' was part of some mesh(es), it will continue to be so,
	// although it is now a triangle
	// but 'new_tri' is not aware of that, for now
	// so we must introduce it to all meshes above 'this'
	new_tri.is_part_of ( *this );

	return new_seg;

} // end of Cell::split in two triangles (cell is a rectangle)


void hidden::for_hanging_nodes_2d ( Mesh & ambient_mesh, Cell & square,
  Cell & AB, Cell * & AB1_p, Cell * & AB2_p, Cell & center, double epsi )

{	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );
	
	// there are three cases
	// 1: the neighbour square is the same size (or even larger)
	//    we need to create a (degenerated) triangle to smoothen the passage
	// 2: the neighbour square has already been split
	//    we destroy the (degenerated) triangle and use the half-segments
	// 3: there is no neighbour square (we are on the boundary of the ambient mesh)
	//    we split the segment but there is no need for a (degenerated) triangle

	enum Situation { create_tri, destroy_tri, boundary };

	Cell & A = AB.base().reverse();
	Cell & B = AB.tip();
	
	Situation situation;
	Cell * cll_p = ambient_mesh.cell_in_front_of ( AB, tag::may_not_exist );
	if ( cll_p == NULL ) situation = boundary;
	else
	{	situation = create_tri;
		if ( ( cll_p->belongs_to(ambient_mesh) ) and
	       ( cll_p->boundary().number_of(tag::cells,tag::of_max_dim) == 3 ) )
		{	// is this the long edge of the triangle ?
			std::map<std::string,void*>::iterator it = cll_p->hook.find("long edge");
			if ( it != cll_p->hook.end() )
			{	Cell * seg_p = (Cell*) (it->second);
				if ( seg_p == AB.hidden_reverse )  situation = destroy_tri;   }          }  }

	switch ( situation ) {
		
	case create_tri : // we need to create a (degenerated) triangle
	{	Cell & middle_of_AB = Cell::point();
		space.interpolate ( middle_of_AB, A, 0.5-epsi, B, 0.5-epsi, center, 2.*epsi );
		Cell & AB1 = Cell::segment ( A.reverse(), middle_of_AB );
		Cell & AB2 = Cell::segment ( middle_of_AB.reverse(), B );
		Cell & AB_tri = Cell::triangle ( AB, AB2.reverse(), AB1.reverse() );
		AB_tri.hook["long edge"] = (void*) (&AB);
		// we must introduce 'AB_tri' to all meshes above 'square'
		AB_tri.is_part_of ( square ); // add_to ?
		AB1_p = & AB1;
		AB2_p = & AB2;                                                                       }
		break;

	case destroy_tri : // we destroy the (degenerated) triangle 
	  // and use the half-segments (yes, this is the long edge of the triangle)
	{	assert ( cll_p != NULL );
		assert ( cll_p->belongs_to (ambient_mesh) );
		Cell & seg1 = cll_p->boundary().cell_in_front_of (A);
		Cell & ver = seg1.tip();
		Cell & seg2 = cll_p->boundary().cell_in_front_of (ver);
		space.interpolate ( ver, A, 0.5, B, 0.5 ); // useless if epsi == 0.
		// when *cll_p is discarded, we do not want seg1 and seg2 to be discarded
		// so we detach them from *cll_p
		seg1.cut_from_bdry_of ( *cll_p );
		seg2.cut_from_bdry_of ( *cll_p );
		// we do not cut AB.reverse() from the boundary of *cll_p
		// thus, when *cll_p is discarded, AB will be discarded, too
		// but we do not want A or B to be discarded so we detach them from AB
		// (alternatively, we could cut AB.reverse() from the boundary of *cll_p
		//  and discard it afterwards; anyway, we have to detach A and B first)
		A.reverse().cut_from_bdry_of ( AB );
		B.cut_from_bdry_of ( AB );
		cll_p->discard(); // removes *cll_p from all meshes above
		AB1_p = & seg1;
		AB2_p = & seg2;                                                        }
		break;

	case boundary : // we are on the boundary of the ambient mesh
	         // we split the segment but there is no need for a new triangle
	{	Cell & AB2 = AB.split ( tag::in_two_segments, tag::cell_is_segment );
		// after split, AB will be the first half of the segment, AB2 the second half
		AB1_p = & AB;
		AB2_p = & AB2;                                                          }
		break;
	} // end of switch

} // end of hidden::for_hanging_nodes_2d


void Cell::split ( const tag::InFourRectangles &, const tag::CellIsRectangle &,
                 const tag::WithIn &, Mesh & ambient_mesh, double epsi         )

// split 'this' rectangle in four smaller rectangles, introducing hanging nodes

// epsi defaults to 0.; represents a slight deformation which turns
// the degenerated triangles visible (less degenerate)

{	// we use the environement manifold
	assert ( Mesh::environment != NULL );
	Manifold & space = * ( Mesh::environment );

	CellIterator it = this->boundary().iter_over ( tag::segments, tag::around );
	it.reset(); assert ( it.in_range() );
	Cell & AB = *it;
	it++; assert ( it.in_range() );
	Cell & BC = *it;
	it++; assert ( it.in_range() );
	Cell & CD = *it;
	it++; assert ( it.in_range() );
	Cell & DA = *it;
	it++; assert ( not it.in_range() );
	Cell & A = DA.tip();
	Cell & B = AB.tip();
	Cell & C = BC.tip();
	Cell & D = CD.tip();

	AB.cut_from_bdry_of (*this);
	BC.cut_from_bdry_of (*this);
	CD.cut_from_bdry_of (*this);
	DA.cut_from_bdry_of (*this);
	// 'this' square has now empty boundary, stays wide open

	Cell & center = Cell::point();
	space.interpolate ( center, A, 0.25, B, 0.25, C, 0.25, D, 0.25 );

	Cell * AB1_p, * AB2_p;
	hidden::for_hanging_nodes_2d ( ambient_mesh, *this, AB, AB1_p, AB2_p, center, epsi );
	Cell & AB1 = * AB1_p;
	Cell & AB2 = * AB2_p;
 	Cell & middle_of_AB = AB1.tip();

	Cell * BC1_p, * BC2_p;
	hidden::for_hanging_nodes_2d ( ambient_mesh, *this, BC, BC1_p, BC2_p, center, epsi );
	Cell & BC1 = * BC1_p;
	Cell & BC2 = * BC2_p;
	Cell & middle_of_BC = BC1.tip();
		
	Cell * CD1_p, * CD2_p;
	hidden::for_hanging_nodes_2d ( ambient_mesh, *this, CD, CD1_p, CD2_p, center, epsi );
	Cell & CD1 = * CD1_p;
	Cell & CD2 = * CD2_p;
	Cell & middle_of_CD = CD1.tip();
		
	Cell * DA1_p, * DA2_p;
	hidden::for_hanging_nodes_2d ( ambient_mesh, *this, DA, DA1_p, DA2_p, center, epsi );
	Cell & DA1 = * DA1_p;
	Cell & DA2 = * DA2_p;
	Cell & middle_of_DA = DA1.tip();
	
	// four new segments :
	Cell & AB_center = Cell::segment ( middle_of_AB.reverse(), center );
	Cell & BC_center = Cell::segment ( middle_of_BC.reverse(), center );
	Cell & CD_center = Cell::segment ( middle_of_CD.reverse(), center );
	Cell & DA_center = Cell::segment ( middle_of_DA.reverse(), center );

	// redefine 'this' rectangle
	AB1.glue_on_bdry_of (*this);
	AB_center.glue_on_bdry_of (*this);
	DA_center.reverse().glue_on_bdry_of (*this);
	DA2.glue_on_bdry_of (*this);

	// three new rectangles
	Cell & sq_B = Cell::rectangle ( AB2, BC1, BC_center, AB_center.reverse() );
	// we must introduce 'sq_B' to all meshes above 'this'
	sq_B.is_part_of (*this);
	Cell & sq_C = Cell::rectangle ( BC2, CD1, CD_center, BC_center.reverse() );
	// we must introduce 'sq_B' to all meshes above 'this'
	sq_C.is_part_of (*this);
	Cell & sq_D = Cell::rectangle ( CD2, DA1, DA_center, CD_center.reverse() );
	// we must introduce 'sq_B' to all meshes above 'this'
	sq_D.is_part_of (*this);

} // end of Cell::split in four rectangles (cell is a rectangle)


