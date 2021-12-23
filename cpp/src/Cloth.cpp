#include "Cloth.h"

namespace csf {

	Cloth::Cloth(
		const std::size_t& nrows_param, 
		const std::size_t& ncols_param, 
		const int& rigidness_param, 
		const double& row_step_param, 
		const double& col_step_param, 
		const double& timestamp_param, 
		const Vector3d& initial_param)
		:	nrows(nrows_param), 
			ncols(ncols_param),
			rigidness(rigidness_param),
			row_step(row_step_param),
			col_step(col_step_param),
			timestamp(timestamp_param),
			initial_position(initial_param)
	{
		particles.reserve(nrows * ncols); // vector for big size, use reserve() first
		double timestamp_squared(timestamp * timestamp); // pre-compute the timestamp_2 then pass it to particle

		/*
		* create particles in a grid, starting from the first particle
		* NB: <= instead of < -- cover all sample points
		* i.e: bmax-bmin = 99.998, bmin=0, bmax=99.998
		* Rounded up: 0 - 100, nrows or ncols = 100
		* therefore if i or j starts from 0, then it can be equal to 100(use "<=")
		*/
		for (std::size_t i = 0; i <= nrows; ++i) {
			for (std::size_t j = 0; j <= ncols; ++j) {
				Vector3d pos(initial_position.v[0] + i * row_step,
							 initial_position.v[1] + j * col_step,
							 initial_position.v[2]);

				particles.emplace_back(Particle(pos, timestamp_squared));
				particles[j + i * ncols].row = i;
				particles[j + i * ncols].col = j;
			}
		}
	}

}
