#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "GroundFilter.h" // Using Point


namespace csf {

	// structure of mypoint
	struct MyPoint {
		double x;
		double y;
		double z;

		MyPoint() :x(0), y(0), z(0) {}
		MyPoint(const double& p_x, const double& p_y, const double& p_z) :
			x(p_x), y(p_y), z(p_z) {}

	};

	// structure of vector
	class Vector3d {
	public:
		double v[3]{}; // the length is known, use the normal array(faster than vector)
	public:
		Vector3d() = default;
		Vector3d(const double& vx, const double& vy, const double& vz) {
			v[0] = vx; v[1] = vy; v[2] = vz;
		}
		Vector3d(const Point& pa, const Point& pb) {
			//v[0] = pb.x - pa.x; 
			//v[1] = pb.y - pa.y; 
			//v[2] = pb.z - pa.z;
			
			v[0] = pb[0] - pa[0]; // x of vector
			v[1] = pb[1] - pa[1]; // y of vector
			v[2] = pb[2] - pa[2]; // z of vector
		} // vector: ab(point to b)

		double vector_length()const;
		Vector3d normalized_vector()const;

		void operator+=(const Vector3d& other);// NB: this operator would change the original vector
		Vector3d operator+(const Vector3d& other)const;
		Vector3d operator-(const Vector3d& other)const;
		Vector3d operator*(const double& s)const;
		Vector3d operator/(const double& s)const;

		double dot_product(const Vector3d& other)const;
		Vector3d cross_product(const Vector3d& other)const;
	};
}

#endif
