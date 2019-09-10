// src/manifem/Util.h 2019.09.08

#ifndef util_h
#define util_h

#include <list>
#include <map>
#include <string>
#include <iostream>

namespace ManiFEM {
	template <class T> class Pointer;
}

// tags are single-object classes (structures, actually).
// Their values are not used; their type is relevant at compile time
// for over-loading functions (mainly factory functions).
// The "regular" convention is that class names and namespaces
// are capitalized. Here we break this convention because we want
// 'tags' to be discrete. For instance, if we wrote 'Tags::reverse_of',
// the reader's eye would catch 'Tags' much before 'reverse_of'.
// Thus, we prefer 'tags::reverse_of'.
// Of course, the user has the choice of 'using namespace tags'.

namespace ManiFEM { namespace tag {
		struct ReverseOf {}; static const ReverseOf reverse_of;
		struct Oriented {}; static const Oriented oriented;
		                    static const Oriented positive;
		struct NotOriented {}; static const NotOriented not_oriented;
		struct Cells {}; static const Cells cells;
		struct Meshes {}; static const Meshes meshes;
		struct Above { }; static const Above above;
		struct OfDim {}; static const OfDim of_dim;
		struct OfMinDim {}; static const OfMinDim of_min_dim;
		struct OfMaxDim {}; static const OfMaxDim of_max_dim;
		struct OfDimPlusOne {}; static const OfDimPlusOne of_dim_plus_one;
		struct OfDimPlusTwo {}; static const OfDimPlusTwo of_dim_plus_two;
		struct OfSameDim {}; static const OfSameDim of_same_dim;
		struct Vertices {}; static const Vertices vertices;
		struct Segments {}; static const Segments segments;
		struct Around {}; static const Around around;
		struct Along {}; static const Along along;
		struct Reverse {}; static const Reverse reverse;
		struct MayNotExist {}; static const MayNotExist may_not_exist;
		struct SurelyExists {}; static const SurelyExists surely_exists;
		struct OnTheFly {}; static const OnTheFly on_the_fly;
		enum WithTriangles { with_triangles, not_with_triangles };
		struct CellIsSegment {}; static const CellIsSegment cell_is_segment;
		struct CellIsRectangle {}; static const CellIsRectangle cell_is_rectangle;
		struct At {}; static const At at;
		struct InTwoSegments {}; static const InTwoSegments in_two_segments;
		struct InTwoTriangles {}; static const InTwoTriangles in_two_triangles;
		struct InFourRectangles {}; static const InFourRectangles in_four_rectangles;
		struct WithIn {}; static const WithIn within;
}	}

using namespace ManiFEM;

template <class T>
class ManiFEM::Pointer

// A tiny raw pointer class, just to avoid conversions to int.

{	public:
	T* p;
	inline Pointer (T* q) : p(q) {};
	inline Pointer (T& q) : p(&q) {};
	inline T& operator*() { return *p; }
	inline T* operator->() { return p; }
};

template < class P, class Q >
inline bool operator== ( Pointer<P> pp, Q q )
{ return pp.p == q; }

template < class P, class Q >
inline bool operator!= ( Pointer<P> pp, Q q )
{ return pp.p != q; }

template < class P, class Q >
inline bool operator== ( Q q, Pointer<P> pp )
{ return q == pp.p; }

template < class P, class Q >
inline bool operator!= ( Q q, Pointer<P> pp )
{ return q != pp.p; }


namespace ManiFEM { namespace hidden {
		template <class T> class VectorOfPointers;
		template <class T> class VectorOfIntegers;
} }

template <class object>
class ManiFEM::hidden::VectorOfPointers

// This is just a tiny wrapper around
// map<short int, list<T>*> and map<short int, map<T,U>*>
// The main contribution is the operator[] which searches for
// a given key and, if it does not exist, creates on-the-fly
// a new list<T> or map<T,U>. And it does this efficiently,
// without unnecessarily repeated searches.
// Inspired in item 24 of the book : Scott Meyers, Effective STL

{ public :
	typedef std::map < short int, object* > index_map;
	typedef typename index_map::iterator iterator;
  typedef typename index_map::reverse_iterator reverse_iterator;
  index_map internal_map;
  object * operator[] ( short int i )
  { //////////////////////////////////////////////////////////////////
    typename index_map::iterator lb = internal_map.lower_bound(i);  //
    if ( (lb==internal_map.end())                                   //
      || (internal_map.key_comp()(i,lb->first)) )                   //
    { // cout << "creating object for key " << i << endl;           //
      typedef typename index_map::value_type VT;                    //
      return internal_map.insert (lb, VT(i,new object))->second; }  //
    else                                                            //
    { // cout << "key " << i << " already exists" << endl;          //
      return lb->second;                                    }       //
//////// code below is functionally equivalent to the above //////////
//  if (internal_map.find(i)==internal_map.end())                   //
//  { typedef typename index_map::value_type VT;                    //
//    return internal_map.insert (VT(i,new object));   }            //
//  else return internal_map[i];                                    //
//////////////////////////////////////////////////////////////////////
  }
  iterator begin () { return internal_map.begin(); }
  iterator end () { return internal_map.end(); }
  reverse_iterator rbegin () { return internal_map.rbegin(); }
  reverse_iterator rend () { return internal_map.rend(); }
};

namespace ManiFEM { namespace hidden {
		template <typename T>
		using IndexedList = hidden::VectorOfPointers < std::list<T> >;
		template <typename T, typename U>
		using IndexedMap = hidden::VectorOfPointers < std::map<T,U> >;
}	}


template <typename object>
class ManiFEM::hidden::VectorOfIntegers

// similarly, a VectorOfIntegers creates on-the-fly a zero value

{ public :
	typedef std::map < short int, object > index_map;
	typedef typename index_map::iterator iterator;
  typedef typename index_map::reverse_iterator reverse_iterator;
  index_map internal_map;
  object & operator[] ( short int i )
  { //////////////////////////////////////////////////////////////////
    typename index_map::iterator lb = internal_map.lower_bound(i);  //
    if ( (lb==internal_map.end())                                   //
      || (internal_map.key_comp()(i,lb->first)) )                   //
    { // cout << "creating object for key " << i << endl;           //
      typedef typename index_map::value_type VT;                    //
			return internal_map.insert(lb,VT(i,(object)(0)))->second;  }  //
    else                                                            //
    { // cout << "key " << i << " already exists" << endl;          //
      return lb->second;                                    }       //
//////// code below is functionally equivalent to the above //////////
//  if (internal_map.find(i)==internal_map.end())                   //
//  { typedef typename index_map::value_type VT;                    //
//    return internal_map.insert (VT(i,new object));   }            //
//  else return internal_map[i];                                    //
//////////////////////////////////////////////////////////////////////
  }
  iterator begin () { return internal_map.begin(); }
  iterator end () { return internal_map.end(); }
  reverse_iterator rbegin () { return internal_map.rbegin(); }
  reverse_iterator rend () { return internal_map.rend(); }
};


/*
int main ()
{
  IndexedList<double> il;
  list<double>* l = il[2];
  list<double>* k = il[2];
  il[1]->push_back(7);
  il[1]->push_back(-7);
  IndexedList<double>::iterator i = il.begin();
  for ( ; i != il.end(); i++)
    cout << i->first << " " << i->second->size() << endl;
  // output is :
  // 1 2
  // 2 0
}
*/

#endif
