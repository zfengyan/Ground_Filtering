/********************************************************************
* GEO1015.2021
* hw03
* Yitong Xia
* 5445825
* Fengyan Zhang
* 5462150
*
* The idea of CSF is from this paper:
* http://www.mdpi.com/2072-4292/8/6/501/htm
* The basic architectures and ideas are inspired by this article (its related code is open source).
*
*********************************************************************/


#include "Particle.h"

namespace csf {

	/*
	* update a particle's position(by gravity)
	* using verlet integration: next_pos = pos * 2 - prev_pos + acc * dt * dt;
	* need to set cur_pos = pre_pos before using(in the constructor of Particle)
	* 
	* The residual of verlet integration is: O(t^4)
	* The precision of timestep: at least 0.01, and the timestep_2 is: 0.0001
	*/
	void Particle::update_position_gravity() {

		if (movable) {
			Vector3d temp(cur_pos);
			cur_pos = cur_pos + cur_pos - pre_pos + acceleration * timestep_2; //update current position
			pre_pos = temp;
		}

	}


	/*
	* update a particle's position(by "virtual spring", obey the Hooke's law)
	* IMPORTANT parameter: rigidness(see "Particle.h")
	* different situations:
	* (1)p1 and p2 all movable: update both
	* (2)p1 movable p2 not: update only p1
	* (3)p2 movable p1 not: update only p2
	* (4)p1 and p2 all unmovable: do nothing
	*/
	void Particle::update_position_spring(const int& rigidness) {
		
		for (std::size_t i = 0; i < neighbours.size(); ++i) {
			Particle* that(neighbours[i]);
			Vector3d height_diff(0, 0, that->cur_pos.v[2] - this->cur_pos.v[2]);

			if (this->movable && that->movable) {

				// set correction vector
				Vector3d correction_vector(height_diff *
					(rigidness > 3 ? 0.5 : doubleMove[rigidness]));

				// adjust the position of the current particle and the neighbour
				this->offset_position(correction_vector);
				that->offset_position(-correction_vector); // update position of p: minus correction_vector
			
			}else if (this->movable && (!that->movable)) {
				Vector3d correction_vector(height_diff *
						(rigidness > 3 ? 1 : singleMove[rigidness]));
				this->offset_position(correction_vector);

			}else if ((!this->movable) && that->movable) {
				Vector3d correction_vector(height_diff *
						(rigidness > 3 ? 1 : singleMove[rigidness]));
				that->offset_position(-correction_vector);

			}

		}
		
	}
}