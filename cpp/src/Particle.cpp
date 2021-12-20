#include "Particle.h"

namespace csf {

	/*
	* update a particle's position(by gravity)
	* using verlet integration:
	* next_pos = pos * 2 - prev_pos + acc * dt * dt;
	*/
	void Particle::update_position_gravity() {
		if (movable) {
			Vector3d temp(cur_pos);
			cur_pos = cur_pos + cur_pos - pre_pos + acceleration * timestamp_2; //update current position
			pre_pos = temp;
		}
	}
}