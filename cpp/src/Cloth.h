#ifndef _CLOTH_H_
#define _CLOTH_H_

/*
* Cloth class
*/

#include "Particle.h"
#include <cmath>

namespace csf {

	class Cloth {
	public:
		std::vector<Particle> particles; // store all particles
		std::size_t nrows; // number of particles in height, x-axis in 3D space
		std::size_t ncols; // number of particles in weight, y-axis in 3D space

		int rigidness; //rigidness of cloth(rigidness: 0,1,2,3 ... )
		double row_step; // depends on the row size of pointcloud and the resolution(user defined)
		double col_step; // depends on the column size of pointcloud and the resolution(user defined)
		double timestep; // each timestep stands for each frame
		
		Vector3d initial_position; // the initial position of the first particle of cloth

	public:

		/*
		* IMPORTANT Construtor 
		*/
		Cloth(
			const std::size_t& nrows_param, 
			const std::size_t& ncols_param,
			const int& rigidness_param, 
			const double& row_step_param,
			const double& col_step_param, 
			const double& timestep_param,
			const Vector3d& initial_param); 

		/*
		* @brief: get the pointer of particle: particels[index]
		* @param: row, col -- the position in the cloth grid
		*/
		Particle* get_particle(const std::size_t& row, const std::size_t& col){
			return &particles[col + row * ncols];
		}


		/*
		* @brief: set the virtual spring of the two particles
		* @param: pointers of two particle
		*/
		void set_virtual_spring(Particle* p1, Particle* p2) {
			p1->neighbours.emplace_back(p2);
			p2->neighbours.emplace_back(p1);
		}


		std::size_t get_cloth_size(){
			return nrows * ncols;
		}


		/*
		* @brief: add force for each particle
		* @param: Vector3d force
		*/
		void addforce_for_particles(const Vector3d& f) {
			for (std::size_t i = 0; i < particles.size(); ++i) {
				particles[i].add_force(f);
			}
		}


		/*
		* @brief:
		* update the status(position) of each particle in the cloth
		* gravity
		*/
		void update_cloth_gravity();


		/*
		* @brief: determine whether the particle hits the terrain
		* compare the Intersection_Height_Value(the height value of corresponding (reversed)LiDAR point)
		* if the current_Height_Value(z coordinate) < Intersect_Height_Value
		* set the z coordinate = Intersect_Height_Value
		* make the particle unmovable
		*/
		void terrain_intersection_check();


		/*
		* @brief:
		* update the status(position) of each particle in the cloth
		* "virtual spring" -- need to execute the intersection_terrain_check() first
		*/
		void update_cloth_spring();


		/*
		* @brief: calculate the maximum height difference
		* between one particle's pre_pos and cur_pos amoung all
		* the cloth particles
		*/
		double calculate_max_diff();


	};

}






#endif // !_CLOTH_H_

