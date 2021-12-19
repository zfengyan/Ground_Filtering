#include "Vector3d.h"
#include <cmath>

namespace csf {
	double Vector3d::vector_length()const {
		return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}

	Vector3d Vector3d::normalized_vector()const {
		const double len(vector_length());
		return len ? Vector3d(v[0] / len, v[1] / len, v[2] / len) : Vector3d();
		// length can't be 0 otherwise return default value
	}

	void Vector3d::operator+=(const Vector3d& other) {
		v[0] += other.v[0];
		v[1] += other.v[1];
		v[2] += other.v[2];
	}

	Vector3d Vector3d::operator+(const Vector3d& other)const {
		return Vector3d(v[0] + other.v[0], v[1] + other.v[1], v[2] + other.v[2]);
	}

	Vector3d Vector3d::operator-(const Vector3d& other) const {
		return Vector3d(v[0] - other.v[0], v[1] - other.v[1], v[2] - other.v[2]);
	}

	Vector3d Vector3d::operator*(const double& s) const {
		return Vector3d(v[0] * s, v[1] * s, v[2] * s);
	}

	Vector3d Vector3d::operator/(const double& s) const {
		return s ? Vector3d(v[0] / s, v[1] / s, v[2] / s) : Vector3d();
	}

	double Vector3d::dot_product(const Vector3d& other) const {
		return (v[0] * other.v[0] + v[1] * other.v[1] + v[2] * other.v[2]);
	}

	Vector3d Vector3d::cross_product(const Vector3d& other) const
	{
		return Vector3d(
			v[1] * other.v[2] - v[2] * other.v[1],
			v[2] * other.v[0] - v[0] * other.v[2],
			v[0] * other.v[1] - v[1] * other.v[0]
		);
	}
}