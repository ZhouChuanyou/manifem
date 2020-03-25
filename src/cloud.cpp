
// example presented in paragraph 8.15 of the manual
// builds a MetricTree over a cloud of points in RR2

#include <iostream>
#include "math.h"
#include "metric-tree.h"

class SqDistanceOnRn

{	public:
	inline double operator() ( const std::vector<double> &, const std::vector<double> & );
};
	
//-----------------------------------------------------------------------------------------------//

inline double SqDistanceOnRn::operator() ( const std::vector<double> & u, const std::vector<double> & v )

{	double res = 0.;
	size_t n = u.size();
	assert ( n == v.size() );
	for ( size_t i = 0; i < n; i++ )
	{	double tmp = u[i] - v[i];
		res += tmp*tmp;           }
	return res;                        }

//-----------------------------------------------------------------------------------------------//


int main ()

{	SqDistanceOnRn sq_dist_Rn;
	MetricTree < std::vector < double >, SqDistanceOnRn > cloud ( sq_dist_Rn, 1., 6. );

	size_t n = 20;
	double step = 2. * 3.14159 / n;
	for ( size_t i = 0; i < n; i++ )
	{	double t = i * step;
		cloud.add ( { std::cos(t), std::sin(t) } );  }
	for ( size_t i = 0; i < n; i++ )
	{	double t = i * step;
		std::list < std::vector < double > > l =
			cloud.find_close_neighbours_of ( { std::cos(t), std::sin(t) }, 0.7 * step * step );
		std::cout << l.size() << std::endl;                                                    }
	// for ( std::list < std::vector < double > > ::const_iterator it = l.begin();
	// 			it != l.end(); it++ )
	//	std::cout << (*it)[0] << " " << (*it)[1] << std::endl;
}
