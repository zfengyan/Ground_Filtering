#include "Particle.h"

namespace csf {

	/*
	* update a particle's position(by gravity)
	* using verlet integration: next_pos = pos * 2 - prev_pos + acc * dt * dt;
	* need to set cur_pos = pre_pos before using(in the constructor of Particle)
	* 
	* The residual of verlet integration is: O(t^4)
	* The precision of timestamp: at least 0.01, and the timestamp_2 is: 0.0001
	*/
	void Particle::update_position_gravity() {

		if (movable) {
			Vector3d temp(cur_pos);
			cur_pos = cur_pos + cur_pos - pre_pos + acceleration * timestamp_2; //update current position
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
		

	}
}