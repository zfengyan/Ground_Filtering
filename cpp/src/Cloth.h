#ifndef _CLOTH_H_
#define _CLOTH_H_

/*
* Cloth class
*/

#include "Particle.h"

namespace csf {

	class Cloth {
	public:
		std::vector<Particle> particles; // store all particles
		std::size_t width; // number of particles in width
		std::size_t height; // number of particles in height

	public:

		/*
		* @brief: get the pointer of particle: particels[index]
		* @param: row, col -- the position in the cloth grid
		*/
		Particle* get_particle(const std::size_t& row, const std::size_t& col)const {
			return &particles[col * width + row];
		}


		/*
		* @brief: set the virtual spring of the two particles
		* @param: pointers of two particle
		*/
		void set_virtual_spring(Particle* p1, Particle* p2) {
			p1->neighbours.emplace_back(p2);
			p2->neighbours.emplace_back(p1);
		}


		std::size_t get_cloth_size()const {
			return width * height;
		}

	};

}






#endif // !_CLOTH_H_

