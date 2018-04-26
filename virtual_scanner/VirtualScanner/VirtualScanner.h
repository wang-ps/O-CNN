#pragma once
#define  SINGLE_PRECISION

#include <vector>
#include <Eigen/Dense>
// CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

using namespace std;
// typedefs
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

typedef Eigen::Matrix<Float, 3, 1>								Vector3F;
typedef Eigen::Matrix<Float, 3, 3>								Matrix3F;
typedef Eigen::Matrix<Float, Eigen::Dynamic, 1> 				VectorXF;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1>  				VectorXb;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1>  					VectorXi;
typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>	MatrixXF;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>		MatrixXi;



typedef CGAL::Simple_cartesian<Float>							K;
typedef K::FT 													FT;
typedef K::Ray_3												Ray;
typedef K::Line_3 												Line;
typedef K::Point_3  											Point;
typedef K::Triangle_3  											Triangle;
typedef std::vector<Triangle>::iterator 						Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator>  			Primitive;
typedef CGAL::AABB_traits<K, Primitive>							AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits>					Tree;
typedef Tree::Point_and_primitive_id							Point_and_primitive_id;
typedef boost::optional< Tree::Intersection_and_primitive_id<Ray>::Type > Ray_intersection;

class VirtualScanner
{
public:
	bool scanning(string filename, int view_num, bool flags, bool normalize);

protected:
	bool save_ply(string filename);
	bool save_binary(string filename, bool flags);
	void calc_views();

protected:
	MatrixXF V_;
	MatrixXi F_;
	MatrixXF N_;
	Vector3F bbMin_, bbMax_;
	vector<Vector3F> view_center_;
	vector<Vector3F> view_dir_;
	vector<Vector3F> dx_;
	vector<Vector3F> dy_;

	vector<float> pts_;
	vector<float> normals_;	
	vector<int> flags_; // indicate whether the normal is reversed

	const int resolution_ = 127;
	const int total_view_num_ = 14;
};

