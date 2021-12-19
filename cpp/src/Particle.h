#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "Vector3d.h"

#define LOWEST_HEIGHT_INF -99999

namespace csf {
	/*
	* Pre-compute the scale factors according to the rigidness(0 to 3)
	* Rigidness: 0,1,2,3

	* For particleMove1, which means one of the two particles is unmovable, then move the other one only:
	* if rigidness = 0: particleMove1 = 0
	* if rigidness = 1: particleMove1 = 0.3, (scale factor of the total distance)
	* if rigidness = 2: particleMove1 = (1-0.3)*0.3+0.3 = 0.51
	* if rigidness = 3: particleMove1 = (1-0.51)*0.3+0.51 = 657

	* For particleMove2, both of the two particles are moved towards each other.
	* if rigidness = 0: particleMove2 = 0
	* if rigidness = 1: particleMove2 = 0.3, (scale factor of the total distance)
	* if rigidness = 2: particleMove2 = (1-0.3*2)*0.3+0.3 = 0.42
	* if rigidness = 3: particleMove2 = (1-0.42*2)*0.3+0.42 = 0.468
	*/
	const double particleMove1[4]{ 0,0.3,0.51,0.657 };
	const double particleMove2[4]{ 0,0.3,0.42,0.468 };

	/*
	* class Particle: a single particle forming the cloth in 3D space
	*/
	class Particle {
	public:
		double mass; // mass of an particle, set to 1
		bool movable; // initialize: true
		Vector3d cur_pos; // current position of the particle
		Vector3d pre_pos; // previous position of the particle, used for iterating positions

		std::size_t row; //position in the cloth grid
		std::size_t col;

		Vector3d acceleration; 
		double timestamp_2; //pre-compute the squared timestamp: timestamp_2

		std::vector<Particle*> neighbours; // neighbour particles of the current particle
		double IntersectionHeightValue; // the height(z value) of the CP(corresponding point in LiDAR point)
		std::size_t correspondLiDAR_index; // the index of corresponding LiDAR point

	public:
		Particle():
			mass(1),
			movable(true),
			cur_pos(Vector3d(0, 0, 0)),
			pre_pos(Vector3d(0, 0, 0)),
			row(0),
			col(0),
			acceleration(Vector3d(0, 0, 0)),
			timestamp_2(0),
			IntersectionHeightValue(LOWEST_HEIGHT_INF),
			correspondLiDAR_index(0){}

		Particle(const Vector3d& pos, const double& timestamp_squared):
			mass(1),
			movable(true),
			cur_pos(pos),
			pre_pos(Vector3d(0, 0, 0)),
			row(0),
			col(0),
			acceleration(Vector3d(0, 0, 0)),
			timestamp_2(timestamp_squared),
			IntersectionHeightValue(LOWEST_HEIGHT_INF),
			correspondLiDAR_index(0){}
	public:

		void add_force(const Vector3d& f) {
			acceleration += f; // a = f/m, mass is assumed as 1
		}

		void offset(const Vector3d& v) {
			if (movable)cur_pos += v;
		}

		void set_unmovable() {
			movable = false;
		}

		/*
		* update a particle's position(by gravity)
		* using verlet integration: 
		* https://www.algorithm-archive.org/contents/verlet_integration/verlet_integration.html
		*/
		void update_position_gravity();

		/*
		* update a particle's position
		* constrained by the "virtual spring" between two particles
		*/
		void update_position_spring();

	};

}

#endif // !_PARTICLE_H_