// manifem/global.cpp 2019.08.17

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "assert.h"
#include "Mesh.h"


Mesh & Mesh::segment_interpolate_yes_or_no  // static
	( Cell & a, Cell & b, size_t n, bool interpolate, NumericField &coord )

{	assert ( a.is_positive() );
	assert ( b.is_positive() );
	Mesh & malha = * ( new Mesh (1) );
	malha.hook["first point"] = (void*) &a;
	malha.hook["last point"] = (void*) &b;
	Cell * prev_point = &a;
	for ( size_t i=1; i < n; ++i )
	{	Cell & p = Cell::point();
		if (interpolate)
		{	double frac = double(i)/double(n);
			coord.interpolate ( p, a, 1.-frac, b, frac );  }
		Cell & seg = Cell::segment ( prev_point->reverse(), p );
		seg.add_to (malha); 
		prev_point = &p;                                         }
	Cell & seg = Cell::segment ( prev_point->reverse(), b );
	seg.add_to (malha);
	return malha;
}

Mesh & Mesh::loop_interpolate_yes_or_no  // static
( size_t n, bool interpolate, NumericField &coord )
// we cannot interpolate actually, how could we ?
{	Mesh & malha = * ( new Mesh (1) );
	Cell & first = Cell::point();
	Cell * prev_point = &first;
	// we prefer do deal always with positive points
	for ( size_t i=1; i < n; ++i )
	{	Cell & p = Cell::point();
		if (interpolate)
		{	double frac = double(i)/double(n);  }
		Cell & seg = Cell::segment ( prev_point->reverse(), p );
		seg.add_to (malha); 
		prev_point = &p;                                         }
	Cell & seg = Cell::segment ( prev_point->reverse(), first );
	seg.add_to (malha);
	return malha;                                                    }


Mesh & Mesh::rectangle_interpolate_yes_or_no  // static
( Mesh & south, Mesh & east, Mesh & north, Mesh & west,
  bool interpolate, NumericField &coord, bool with_triangles )

{ // recover the corners from the sides
	Cell *SE, *NE, *SW, *NW;
	size_t N_horiz, N_vert;
	// we get the first and last points of each side using freshly reset iterators
	// see paragraph 5.2 in the manual
	{ // just to keep it1 and it2 local
	CellIterator it1 = south.iter_over ( tag::vertices, tag::along );
	it1.reset(); SW = &(*it1);
	CellIterator it2 = south.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); SE = &(*it2);
	if ( south.is_positive() )
		N_horiz = south.cells[1]->size();
	else
		N_horiz = south.reverse().cells[1]->size();
	} { // just to keep it1 and it2 local
	CellIterator it1 = north.iter_over ( tag::vertices, tag::along);
	it1.reset(); NE = &(*it1);
	CellIterator it2 = north.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); NW = &(*it2);
	if ( north.is_positive() )
		assert ( N_horiz == north.cells[1]->size() );
	else
		assert ( N_horiz == north.reverse().cells[1]->size() );
	} { // just to keep it1 and it2 local
	CellIterator it1 = east.iter_over ( tag::vertices, tag::along );
	it1.reset(); assert ( SE == &(*it1) );
	CellIterator it2 = east.iter_over ( tag::vertices, tag::along, tag::reverse );
	it2.reset(); assert ( NE == &(*it2) );
	if ( east.is_positive() )
		N_vert = east.cells[1]->size();
	else
		N_vert = east.reverse().cells[1]->size();
	} { // just to keep it1 and it2 local
	CellIterator it1 = west.iter_over ( tag::vertices, tag::along );
	it1.reset(); assert ( NW == &(*it1) );
	CellIterator it2 = west.iter_over ( tag::vertices, tag::along, tag::reverse);
	it2.reset(); assert ( SW == &(*it2) );
	if ( west.is_positive() )
		assert ( N_vert == west.cells[1]->size() );
	else
		assert ( N_vert == west.reverse().cells[1]->size() );
	} // just to keep it1 and it2 local
		
	// prepare the mesh
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
	
	// prepare ground
	std::list <Cell*> ground ;
	Cell * P = SW;
	for ( size_t j = 1; j <= N_horiz; ++j )
	{	Cell * s = south.cell_in_front_of ( P, tag::may_not_exist );
		assert ( s != NULL );
		ground.push_back (s);
		P = & ( s->tip() );                                              }

	// start mesh generation
	for ( size_t i = 1; i < N_vert; ++i )
	{	std::list <Cell*>::iterator it = ground.begin();
		Cell & A = (*it)->base().reverse();
		Cell * s4_p = west.cell_behind ( &A, tag::may_not_exist );
		assert ( s4_p != NULL );
		Cell * D_p = & (  s4_p->base().reverse() );
		for ( size_t j = 1; j < N_horiz; j++ )
		{	//Cell & s1 = south.cell_in_front_of (A);
			Cell * s1_p = *it;
			Cell & B = s1_p->tip();
			Cell & C = Cell::point();
			if (interpolate)
			{	double frac_N = double(i) / double(N_vert),
					frac_E = double(j) / double(N_horiz),
					frac_SW = (1.-frac_N) * (1.-frac_E),
					frac_SE = (1.-frac_N) * frac_E,
					frac_NE = frac_N * frac_E,
					frac_NW = frac_N * (1.-frac_E);
				coord.interpolate ( C, *NE, frac_NE, *NW, frac_NW, *SE, frac_SE, *SW, frac_SW);  }
			Cell & s2 = Cell::segment ( B.reverse(), C );
			Cell & s3 = Cell::segment ( C.reverse(), *D_p );
			if ( with_triangles )
			{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
				Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
				Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
				T1.add_to (rectangle);
				T2.add_to (rectangle);                                 }
			else // with quadrilaterals
			{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p );
				Q.add_to (rectangle);                                 }
			// 's3' is on the ceiling; we keep it in the 'ground' list, since
			// it will be on the ground when we build the next layer of cells
			*it = & ( s3.reverse() ); // 'it' points into the 'ground' list
			it++;
			D_p = &C;
			s4_p = & ( s2.reverse() );
		} // end of if (interpolate)
		//Cell & s1 = south.cell_in_front_of (A);
		Cell * s1_p = *it;
		Cell & B = s1_p->tip();
		Cell & s2 = east.cell_in_front_of (B);
		Cell & C = s2.tip();
		Cell & s3 = Cell::segment ( C.reverse(), *D_p );
		if ( with_triangles )
		{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
			Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
			Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
			T1.add_to (rectangle);
			T2.add_to (rectangle);                                 }
		else // with quadrilaterals
		{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p );
			Q.add_to (rectangle);                                 }
		*it = & ( s3.reverse() );
	} // end of for j
	std::list <Cell*>::iterator it = ground.begin();
	Cell * s4_p = west.cell_in_front_of ( NW, tag::may_not_exist );
	assert ( s4_p != NULL );
	Cell * D_p = NW;
	for (size_t j=1; j < N_horiz; j++)
	{	//Cell & s1 = south.cell_in_front_of (A);
		Cell * s1_p = *it;
		Cell & B = s1_p->tip();
		Cell & s3 = north.cell_behind ( *D_p );
		Cell & C = s3.base().reverse();
		Cell & s2 = Cell::segment ( B.reverse(), C );
		if ( with_triangles )
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
	//Cell & s1 = south.cell_in_front_of (A);
	Cell * s1_p = *it;
	Cell & B = s1_p->tip();
	Cell & s2 = east.cell_in_front_of (B);
	Cell & C = s2.tip();
	assert ( &C == NE );
	Cell & s3 = north.cell_behind ( *D_p );
	if ( with_triangles )
	{	Cell & BD = Cell::segment ( B.reverse(), *D_p );
		Cell & T1 = Cell::triangle ( BD.reverse(), s2, s3 );
		Cell & T2 = Cell::triangle ( BD, *s4_p, *s1_p );
		T1.add_to (rectangle);
		T2.add_to (rectangle);                                 }
	else // with quadrilaterals
	{	Cell & Q = Cell::rectangle ( *s1_p, s2, s3, *s4_p);
		Q.add_to (rectangle);                                 }
	return rectangle;

} // end of Mesh::rectangle_interpolate_yes_or_no


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
		CellIterator itt = m.iter_over ( tag::cells_of_max_dim, tag::oriented );
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
	{	CellIterator it1 = cell1_pos->boundary().iter_over ( tag::cells_of_max_dim, tag::oriented );
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
		CellIterator it2 = cell2_pos->boundary().iter_over ( tag::cells_of_max_dim, tag::oriented );
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

	CellIterator it1 = mesh1.iter_over ( tag::cells_of_max_dim, tag::oriented );	
	CellIterator it2 = mesh2.iter_over ( tag::cells_of_max_dim, tag::oriented );

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

{	assert ( Mesh::environment != NULL );
	NumericField &coord = * ( Mesh::environment->coord_field );

	assert (coord.size() == 2);
	// we split 'coord' into its components
	OneDimField &x = coord[0], &y = coord[1];
	
	double xmin, xmax, ymin, ymax, maxside;
	
	{ // just a block for hiding variables
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0, tag::not_oriented );
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
	file_ps << "%%BoundingBox:  " << -border*scale_factor << " " << -border*scale_factor << " "
					<< scale_factor*(xmax+border) + translation_x << "   "
					<< scale_factor*(ymax+border) + translation_y << std::endl;
	file_ps << "%%EndComments" << std::endl << std::endl;
	
	file_ps << "gsave" << std::endl;
	// file_ps << "/m{moveto}def" << std::endl;
	// file_ps << "/l{lineto}def" << std::endl;
	// file_ps << "/s{stroke}def" << std::endl;
	file_ps << translation_x << " " << translation_y << " translate" << std::endl;
	file_ps << scale_factor << " " << scale_factor << " scale" << std::endl << std::endl;

	file_ps << "gsave 0.003 setlinewidth" << std::endl;
	
	{ // just a block for hiding variables
	CellIterator it = this->iter_over ( tag::cells_of_dim, 1, tag::not_oriented );
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
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0, tag::not_oriented );
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


void Mesh::export_msh ( std::string f )

// the numbering of vertices is produced on-the-fly

{	assert ( Mesh::environment != NULL );
	NumericField &coord = * ( Mesh::environment->coord_field );

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl;
	
  file_msh << "$Nodes" << std::endl << this->cells[0]->size() << std::endl;

	std::map < Cell*, size_t > numbering;

	{ // just to make variables local : it, counter, x, y
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0, tag::not_oriented );
	size_t counter = 0;
	OneDimField & x = coord[0], & y = coord[1];
	if (coord.size() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	++counter;
			Cell & p = *it;
			numbering[&p] = counter;
			file_msh << counter << " " << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	if ( coord.size() == 3 )
	{	OneDimField &z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	++counter;
			Cell & p = *it;
			numbering[&p] = counter;
			file_msh << counter << " " << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl;
	} // just to make variables local : it, counter, x, y

	file_msh << "$Elements" << std::endl;
	file_msh << this->number_of ( tag::cells_of_max_dim ) << std::endl;

	if ( this->dim == 2 )
	{	CellIterator it = this->iter_over ( tag::cells_of_max_dim, tag::oriented );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell & elem = *it;
			if ( elem.boundary().number_of ( tag::cells_of_max_dim ) == 3 ) // a triangle
				file_msh << counter << " 2 0 ";
			else // a quadrilateral
			{	assert ( elem.boundary().number_of ( tag::cells_of_max_dim ) == 4 );
				file_msh << counter << " 3 0 ";                            }
			CellIterator itt = elem.boundary().iter_over ( tag::vertices, tag::around );
			for ( itt.reset(); itt.in_range(); ++itt )
			{	Cell & p = *itt;
				file_msh << numbering[&p] << " ";   }
			file_msh << std::endl;                                                            }  }
	else
	{	assert ( this->dim == 3);
		CellIterator it = this->iter_over ( tag::cells_of_max_dim, tag::oriented );
		size_t counter = 0;
		for ( it.reset() ; it.in_range(); it++)
		{	++counter;
			Cell & elem = *it;
			size_t n_faces = elem.boundary().number_of ( tag::cells_of_max_dim );
			if ( n_faces == 4 ) // a tetrahedron
				file_msh << counter << " 4 0 ";  // to finish !
			else if ( n_faces == 6 )
			{	// 3d parallelogram = 8-node hexaedron = cube
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 5 0 ";
				CellIterator itt = elem.boundary().iter_over ( tag::cells_of_max_dim, tag::oriented );
				itt.reset();
				Cell & back = *itt; // square face behind the cube
				// back is 0321 in gmsh's documentation
				assert ( back.boundary().number_of ( tag::cells_of_max_dim ) == 4 );
				CellIterator itv = back.boundary().iter_over
					( tag::vertices, tag::around, tag::reverse );
				// reverse because we want the vertices ordered as 0, 1, 2, 3
				itv.reset();  Cell & ver_0 = *itv;
				Cell & seg_03 = back.boundary().cell_in_front_of(ver_0);
				for ( ; itv.in_range(); ++itv )
				{	Cell & p = *itv;
					file_msh << numbering[&p] << " ";   }
				Cell & left_wall = elem.boundary().cell_in_front_of(seg_03); // square face on the left
				// left_wall is 0473 in gmsh's documentation
				assert ( left_wall.boundary().number_of ( tag::cells_of_max_dim ) == 4 );
				Cell & seg_04 = left_wall.boundary().cell_in_front_of(ver_0);
				Cell & ver_4 = seg_04.tip();
				Cell & seg_47 = left_wall.boundary().cell_in_front_of(ver_4);
				Cell & front = elem.boundary().cell_in_front_of(seg_47); // square face in front
				// front is 4567 in gmsh's documentation
				assert ( front.boundary().number_of ( tag::cells_of_max_dim ) == 4 );
				CellIterator itvv = front.boundary().iter_over ( tag::vertices, tag::around );
				itvv.reset(ver_4);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell & p = *itvv;
					file_msh << numbering[&p] << " ";   }                                                  }
			else
			{	assert( n_faces == 5 );
				// triangular prism = 6-node prism
				// see http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
				// and http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
				file_msh << counter << " 6 0 ";
				CellIterator itt = elem.boundary().iter_over ( tag::cells_of_max_dim, tag::oriented );
				short int n_tri = 0, n_rect = 0;
				Cell* base_p;
				for( itt.reset(); itt.in_range(); ++itt )
				{	Cell & face = *itt; // every face elem
			    size_t n_edges = face.boundary().number_of ( tag::cells_of_max_dim );
					if ( n_edges == 3 )  { n_tri++;   base_p = &face;           }
					else                 { n_rect++;  assert ( n_edges == 4 );  }   }
				assert ( n_tri == 2 );  assert ( n_rect == 3 );
				// base is 021 in gmsh's documentation
				assert ( base_p->boundary().number_of ( tag::cells_of_max_dim ) == 3 );
				CellIterator itv = base_p->boundary().iter_over
					( tag::vertices, tag::around, tag::reverse );
				// we set it reverse because we want the vertices ordered as 0, 1, 2
				itv.reset();
				Cell & ver_0 = *itv;
				Cell & seg_02 = base_p->boundary().cell_in_front_of(ver_0);
				for (  ; itv.in_range(); ++itv )
				{	Cell & p = *itv;
					file_msh << numbering[&p] << " ";   }
				Cell & right_wall = elem.boundary().cell_in_front_of(seg_02);
				// right_wall is 0352 in gmsh's documentation
				assert ( right_wall.boundary().number_of ( tag::cells_of_max_dim ) == 4 );
				Cell & seg_03 = right_wall.boundary().cell_in_front_of(ver_0);
				Cell & ver_3 = seg_03.tip();
				Cell & seg_35 = right_wall.boundary().cell_in_front_of(ver_3);
				Cell & roof = elem.boundary().cell_in_front_of(seg_35); //
				// roof is 345 in gmsh's documentation
				assert ( roof.boundary().number_of ( tag::cells_of_max_dim ) == 3 );
				CellIterator itvv = roof.boundary().iter_over ( tag::vertices, tag::around );
				itvv.reset(ver_3);
				for ( ; itvv.in_range(); ++itvv )
				{	Cell & p = *itvv;
					file_msh << numbering[&p] << " ";   }                                                }
			file_msh << std::endl;                                                                      } }
	file_msh << "$EndElements" << std::endl;
	
	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


void Mesh::export_msh ( std::string f, FiniteElement & fe )

// the numbering of vertices is given by the degrees of freedom provided by 'fe'

{	assert ( Mesh::environment != NULL );
	NumericField &coord = * ( Mesh::environment->coord_field );

	std::ofstream file_msh (f);
	file_msh << "$MeshFormat" << std::endl << "2.2 0 8" << std::endl;
	file_msh << "$EndMeshFormat" << std::endl << "$Nodes" << std::endl;
	file_msh << this->cells[0]->size() << std::endl;

	std::map < Cell*, int > numbering;

	{ // just to make variables local : it, x, y
	CellIterator it = this->iter_over ( tag::cells_of_dim, 0, tag::not_oriented );
	OneDimField & x = coord[0], & y = coord[1];
	if (coord.size() == 2)
	{	for ( it.reset() ; it.in_range(); it++ )
		{	Cell & p = *it;
			file_msh << fe.get_dof(p)+1 << " " << x(p) << " " << y(p) << " " << 0 << std::endl;  }  }
	if (coord.size() == 3)
	{	OneDimField &z = coord[2];
		for ( it.reset() ; it.in_range(); it++ )
		{	Cell & p = *it;
			file_msh << fe.get_dof(p)+1 << " " << x(p) << " " << y(p) << " " << z(p) << std::endl;  }  }
	file_msh << "$EndNodes" << std::endl << "$Elements" << std::endl;
	file_msh << this->cells[2]->size() << std::endl;
	} // just to make variables local : it, x, y

	{ // just to make variables local : it, counter
	CellIterator it = this->iter_over ( tag::cells_of_max_dim, tag::oriented );
	size_t counter = 0;
	for ( it.reset() ; it.in_range(); it++)
	{	++counter;
		Cell & elem = *it;
		if ( elem.boundary().cells[1]->size() == 3 )
			file_msh << counter << " 2 0 ";
		else
		{	assert ( elem.boundary().cells[1]->size() == 4 );
			file_msh << counter << " 3 0 ";                    }
		CellIterator itt = elem.boundary().iter_over ( tag::vertices, tag::around );
		itt.reset();
		size_t s = elem.boundary().cells[0]->size();
		for ( size_t i = 0; i < s; ++i )
		{	assert ( itt.in_range() );
			Cell & p = *itt;
			file_msh << fe.get_dof(p)+1 << " ";
			itt++;                               }
		assert ( ! ( itt.in_range() ) );
		file_msh << std::endl;                                                }
	file_msh << "$EndElements" << std::endl;
	} // just to make variables local : it, counter

	if ( ! file_msh.good() )
	{	std::cerr << "error writing msh file" << std::endl;
		exit (1);                                    }

} // end of Mesh::export_msh


Cell & Cell::cut ( const tag::At &, Cell & ver, const tag::CellIsARectangle & )

// cut 'this' rectangle along a diagonal
// which diagonal ? the one passing through vertex 'ver'
// returns the newly created segment, beginning at 'ver'

// we must also cut the reverse cell, if it exists - to do

{
	Mesh & bdry = this->boundary();
	// let's make sure 'this' is a rectangle (quadrilateral, actually)
	assert ( this->dim == 2 );
	assert ( bdry.cells[0]->size() == 4 );
	assert ( bdry.cells[1]->size() == 4 );

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
	MeshIterator it = this->iter_over ( tag::meshes_above, tag::oriented );
	for ( it.reset(); it.in_range(); it++ ) new_tri.add_to (*it);

	return new_seg;
}


void Cell::cut ( const tag::InFour &, const tag::CellIsARectangle &,
                 const tag::WithIn &, Mesh & ambient_mesh, double epsi )

// cuts a rectangle in four smaller rectangles, introducing hanging nodes

// epsi defaults to 0.; represents a slight deformation which turns
// the degenerated triangles visible

{	NumericField * coord = Mesh::environment->coord_field;

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
	coord->interpolate ( center, A, 0.25, B, 0.25, C, 0.25, D, 0.25 );

	// there are three cases
	// 1: the neighbour square is the same size (or even larger)
	//    we need to create a (degenerated) triangle to smoothen the passage
	// 2: the neighbour square has already been split
	//    we destroy the (degenerated) triangle and use the half-segments
	// 3: there is no neighbour square
	//    we split the segment but there is no need for a (degenerated) triangle

	Cell * middle_of_AB_p, * AB1_p, * AB2_p;
	short int situation;
	Cell * cll = ambient_mesh.cell_in_front_of ( &AB, tag::may_not_exist );
	if ( cll == NULL ) situation = 3;
	else
	{	situation = 1;
		if ( ( cll->belongs_to ( ambient_mesh ) ) and
	       ( cll->boundary().cells[1]->size() == 3 ) )
		{	// is this the long edge of the triangle ?
			std::map<std::string,void*>::iterator it = cll->hook.find("long edge");
			if ( it != cll->hook.end() )
			{	Cell * seg_p = (Cell*) (it->second);
				if ( seg_p == AB.hidden_reverse )
					situation = 2;                      }                                }  }
	switch ( situation ) {
	case 1 : // we need to create a (degenerated) triangle
	{	Cell & middle_of_AB = Cell::point();
		coord->interpolate ( middle_of_AB, A, 0.5-epsi, B, 0.5-epsi, center, 2.*epsi );
		Cell & AB1 = Cell::segment ( A.reverse(), middle_of_AB );
		Cell & AB2 = Cell::segment ( middle_of_AB.reverse(), B );
		Cell & AB_tri = Cell::triangle ( AB, AB2.reverse(), AB1.reverse() );
		AB_tri.hook["long edge"] = (void*) (&AB);
		std::cout << "in Cell::cut, adding cell with " << AB_tri.boundary().cells[1]->size()
							<< " sides" << std::endl;
		AB_tri.add_to ( ambient_mesh );
		middle_of_AB_p = & middle_of_AB;
		AB1_p = & AB1;
		AB2_p = & AB2;                                                         }
		break;
	case 2 : // we destroy the (degenerated) triangle and use the half-segments
		// is this the long edge of the triangle ?
	{	assert ( cll != NULL );
		assert ( cll->belongs_to (ambient_mesh) );
		Cell & seg1 = cll->boundary().cell_in_front_of (A);
		Cell & ver = seg1.tip();
		Cell & seg2 = cll->boundary().cell_in_front_of (ver);
		coord->interpolate ( ver, A, 0.5, B, 0.5 );
		std::cout << "destroying cell with " << cll->boundary().cells[1]->size()
							<< " sides, at AB" << std::endl;
		AB.reverse().cut_from_bdry_of ( *cll );
		seg1.cut_from_bdry_of ( *cll );
		seg2.cut_from_bdry_of ( *cll );
		cll->remove_from ( ambient_mesh );
		cll->discard();
		A.reverse().cut_from_bdry_of ( AB );
		B.cut_from_bdry_of ( AB );
		AB.discard();
		middle_of_AB_p = & ver;
		AB1_p = & seg1;
		AB2_p = & seg2;                                                        }
		break;
	case 3 : // we split the segment but there is no need for a new triangle
	{	Cell & middle_of_AB = Cell::point();
		coord->interpolate ( middle_of_AB, A, 0.5, B, 0.5 );
		B. cut_from_bdry_of ( AB );
		middle_of_AB. glue_on_bdry_of ( AB );
		Cell & AB2 = Cell::segment ( middle_of_AB.reverse(), B );
		middle_of_AB_p = & middle_of_AB;
		AB1_p = & AB;
		AB2_p = & AB2;                                                          }
		break;
	} // end of switch
	Cell & middle_of_AB = * middle_of_AB_p;
	Cell & AB1 = * AB1_p;
	Cell & AB2 = * AB2_p;

	Cell * middle_of_BC_p, * BC1_p, * BC2_p;
	cll = ambient_mesh.cell_in_front_of ( &BC, tag::may_not_exist );
	if ( cll == NULL ) situation = 3;
	else
	{	situation = 1;
		if ( ( cll->belongs_to ( ambient_mesh ) ) and
	       ( cll->boundary().cells[1]->size() == 3 ) )
		{	// is this the long edge of the triangle ?
			std::map<std::string,void*>::iterator it = cll->hook.find("long edge");
			if ( it != cll->hook.end() )
			{	Cell * seg_p = (Cell*) (it->second);
				if ( seg_p == BC.hidden_reverse )
					situation = 2;                      }                                  }  }
	switch ( situation ) {
	case 1 : // we need to create a (degenerated) triangle
	{	Cell & middle_of_BC = Cell::point();
		coord->interpolate ( middle_of_BC, B, 0.5-epsi, C, 0.5-epsi, center, 2.*epsi );
		Cell & BC1 = Cell::segment ( B.reverse(), middle_of_BC );
		Cell & BC2 = Cell::segment ( middle_of_BC.reverse(), C );
		Cell & BC_tri = Cell::triangle ( BC, BC2.reverse(), BC1.reverse() );
		BC_tri.hook["long edge"] = (void*) (&BC);
		std::cout << "in Cell::cut, adding cell with " << BC_tri.boundary().cells[1]->size()
							<< " sides" << std::endl;
		BC_tri.add_to ( ambient_mesh );
		middle_of_BC_p = & middle_of_BC;
		BC1_p = & BC1;
		BC2_p = & BC2;                                                         }
		break;
	case 2 : // we destroy the (degenerated) triangle and use the half-segments
		// is this the long edge of the triangle ?
	{	assert ( cll != NULL );
		assert ( cll->belongs_to (ambient_mesh) );
		Cell & seg1 = cll->boundary().cell_in_front_of (B);
		Cell & ver = seg1.tip();
		Cell & seg2 = cll->boundary().cell_in_front_of (ver);
		coord->interpolate ( ver, B, 0.5, C, 0.5 );
		std::cout << "destroying cell with " << cll->boundary().cells[1]->size()
							<< " sides, at BC" << std::endl;
		BC.reverse().cut_from_bdry_of ( *cll );
		seg1.cut_from_bdry_of ( *cll );
		seg2.cut_from_bdry_of ( *cll );
		cll->remove_from ( ambient_mesh );
		cll->discard();
		B.reverse().cut_from_bdry_of ( BC );
		C.cut_from_bdry_of ( BC );
		BC.discard();
		middle_of_BC_p = & ver;
		BC1_p = & seg1;
		BC2_p = & seg2;                                                        }
		break;
	case 3 : // we split the segment but there is no need for a new triangle
	{	Cell & middle_of_BC = Cell::point();
		coord->interpolate ( middle_of_BC, B, 0.5, C, 0.5 );
		C. cut_from_bdry_of ( BC );
		middle_of_BC. glue_on_bdry_of ( BC );
		Cell & BC2 = Cell::segment ( middle_of_BC.reverse(), C );
		middle_of_BC_p = & middle_of_BC;
		BC1_p = & BC;
		BC2_p = & BC2;                                                         }
		break;
	} // end of switch
	Cell & middle_of_BC = * middle_of_BC_p;
	Cell & BC1 = * BC1_p;
	Cell & BC2 = * BC2_p;
		
	Cell * middle_of_CD_p, * CD1_p, * CD2_p;
	cll = ambient_mesh.cell_in_front_of ( &CD, tag::may_not_exist );
	if ( cll == NULL ) situation = 3;
	else
	{	situation = 1;
		if ( ( cll->belongs_to ( ambient_mesh ) ) and
	       ( cll->boundary().cells[1]->size() == 3 ) )
		{	// is this the long edge of the triangle ?
			std::map<std::string,void*>::iterator it = cll->hook.find("long edge");
			if ( it != cll->hook.end() )
			{	Cell * seg_p = (Cell*) (it->second);
				if ( seg_p == CD.hidden_reverse )
					situation = 2;                      }                                }  }
	switch ( situation ) {
	case 1 : // we need to create a (degenerated) triangle
	{	Cell & middle_of_CD = Cell::point();
		coord->interpolate ( middle_of_CD, C, 0.5-epsi, D, 0.5-epsi, center, 2.*epsi );
		Cell & CD1 = Cell::segment ( C.reverse(), middle_of_CD );
		Cell & CD2 = Cell::segment ( middle_of_CD.reverse(), D );
		Cell & CD_tri = Cell::triangle ( CD, CD2.reverse(), CD1.reverse() );
		CD_tri.hook["long edge"] = (void*) (&CD);
		std::cout << "in Cell::cut, adding cell with " << CD_tri.boundary().cells[1]->size()
							<< " sides" << std::endl;
		CD_tri.add_to ( ambient_mesh );
		middle_of_CD_p = & middle_of_CD;
		CD1_p = & CD1;
		CD2_p = & CD2;                                                         }
		break;
	case 2 : // we destroy the (degenerated) triangle and use the half-segments
		// is this the long edge of the triangle ?
	{	assert ( cll != NULL );
		assert ( cll->belongs_to (ambient_mesh) );
		Cell & seg1 = cll->boundary().cell_in_front_of (C);
		Cell & ver = seg1.tip();
		Cell & seg2 = cll->boundary().cell_in_front_of (ver);
		coord->interpolate ( ver, C, 0.5, D, 0.5 );
		std::cout << "destroying cell with " << cll->boundary().cells[1]->size()
							<< " sides, at CD" << std::endl;
		CD.reverse().cut_from_bdry_of ( *cll );
		seg1.cut_from_bdry_of ( *cll );
		seg2.cut_from_bdry_of ( *cll );
		cll->remove_from ( ambient_mesh );
		cll->discard();
		std::cout << "aha" << std::endl << std::flush;
		C.reverse().cut_from_bdry_of ( CD );
		D.cut_from_bdry_of ( CD );
		CD.discard();
		middle_of_CD_p = & ver;
		CD1_p = & seg1;
		CD2_p = & seg2;                                                         }
		break;
	case 3 : // we split the segment but there is no need for a new triangle
	{	Cell & middle_of_CD = Cell::point();
		coord->interpolate ( middle_of_CD, C, 0.5, D, 0.5 );
		D. cut_from_bdry_of ( CD );
		middle_of_CD. glue_on_bdry_of ( CD );
		Cell & CD2 = Cell::segment ( middle_of_CD.reverse(), D );
		middle_of_CD_p = & middle_of_CD;
		CD1_p = & CD;
		CD2_p = & CD2;                                                         }
		break;
	} // end of switch
	Cell & middle_of_CD = * middle_of_CD_p;
	Cell & CD1 = * CD1_p;
	Cell & CD2 = * CD2_p;
		
	Cell * middle_of_DA_p, * DA1_p, * DA2_p;
	cll = ambient_mesh.cell_in_front_of ( &DA, tag::may_not_exist );
	if ( cll == NULL ) situation = 3;
	else
	{	situation = 1;
		if ( ( cll->belongs_to ( ambient_mesh ) ) and
	       ( cll->boundary().cells[1]->size() == 3 ) )
		{	// is this the long edge of the triangle ?
			std::map<std::string,void*>::iterator it = cll->hook.find("long edge");
			if ( it != cll->hook.end() )
			{	Cell * seg_p = (Cell*) (it->second);
				if ( seg_p == DA.hidden_reverse )
					situation = 2;                      }                               }  }
	switch ( situation ) {
	case 1 : // we need to create a (degenerated) triangle
	{	Cell & middle_of_DA = Cell::point();
		coord->interpolate ( middle_of_DA, D, 0.5-epsi, A, 0.5-epsi, center, 2.*epsi );
		Cell & DA1 = Cell::segment ( D.reverse(), middle_of_DA );
		Cell & DA2 = Cell::segment ( middle_of_DA.reverse(), A );
		Cell & DA_tri = Cell::triangle ( DA, DA2.reverse(), DA1.reverse() );
		DA_tri.hook["long edge"] = (void*) (&DA);
		std::cout << "in Cell::cut, adding cell with " << DA_tri.boundary().cells[1]->size()
							<< " sides" << std::endl;
		DA_tri.add_to ( ambient_mesh );
		middle_of_DA_p = & middle_of_DA;
		DA1_p = & DA1;
		DA2_p = & DA2;                                                         }
		break;
	case 2 : // we destroy the (degenerated) triangle and use the half-segments
		// is this the long edge of the triangle ?
	{	assert ( cll != NULL );
		assert ( cll->belongs_to (ambient_mesh) );
		Cell & seg1 = cll->boundary().cell_in_front_of (D);
		Cell & ver = seg1.tip();
		Cell & seg2 = cll->boundary().cell_in_front_of (ver);
		coord->interpolate ( ver, D, 0.5, A, 0.5 );
		std::cout << "destroying cell with " << cll->boundary().cells[1]->size()
							<< " sides, at DA" << std::endl;
		DA.reverse().cut_from_bdry_of ( *cll );
		seg1.cut_from_bdry_of ( *cll );
		seg2.cut_from_bdry_of ( *cll );
		cll->remove_from ( ambient_mesh );
		cll->discard();
		D.reverse().cut_from_bdry_of ( DA );
		A.cut_from_bdry_of ( DA );
		DA.discard();
		middle_of_DA_p = & ver;
		DA1_p = & seg1;
		DA2_p = & seg2;                                                        }
		break;
	case 3 : // we split the segment but there is no need for a new triangle
	{	Cell & middle_of_DA = Cell::point();
		coord->interpolate ( middle_of_DA, D, 0.5, A, 0.5 );
		A. cut_from_bdry_of ( DA );
		middle_of_DA. glue_on_bdry_of ( DA );
		Cell & DA2 = Cell::segment ( middle_of_DA.reverse(), A );
		middle_of_DA_p = & middle_of_DA;
		DA1_p = & DA;
		DA2_p = & DA2;                                                         }
		break;
	} // end of switch
	Cell & middle_of_DA = * middle_of_DA_p;
	Cell & DA1 = * DA1_p;
	Cell & DA2 = * DA2_p;
	
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
	std::cout << "in Cell::cut, keeping cell with " << this->boundary().cells[1]->size()
						<< " sides" << std::endl;

	// three new rectangles
	Cell & sq_B = Cell::rectangle ( AB2, BC1, BC_center, AB_center.reverse() );
	std::cout << "in Cell::cut, adding cell with " << sq_B.boundary().cells[1]->size()
						<< " sides" << std::endl;
	sq_B.add_to ( ambient_mesh );
	Cell & sq_C = Cell::rectangle ( BC2, CD1, CD_center, BC_center.reverse() );
	std::cout << "in Cell::cut, adding cell with " << sq_C.boundary().cells[1]->size()
						<< " sides" << std::endl;
	sq_C.add_to ( ambient_mesh );
	Cell & sq_D = Cell::rectangle ( CD2, DA1, DA_center, CD_center.reverse() );
	std::cout << "in Cell::cut, adding cell with " << sq_D.boundary().cells[1]->size()
						<< " sides" << std::endl;
	sq_D.add_to ( ambient_mesh );

} // end of Cell::cut

