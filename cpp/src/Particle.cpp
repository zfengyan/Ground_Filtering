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
}