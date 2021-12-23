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
		* Rounded up: 0 - 100, nrows or ncols = 100 + 1 = 101
		*/

		for (std::size_t i = 0; i < nrows; ++i) {
			for (std::size_t j = 0; j < ncols; ++j) {
				Vector3d pos(initial_position.v[0] + i * row_step,
							 initial_position.v[1] + j * col_step,
							 initial_position.v[2]);

				particles.emplace_back(Particle(pos, timestamp_squared));
				particles[j + i * ncols].row = i;
				particles[j + i * ncols].col = j;
			}
		}


		/*
		* Set the virtual spring connection between (directly)adjacent neighbors
		* Define 'block': every small square formed with 4 particles
		* For each 'block', starting from the up-left corner particle: (row, col)
		* Connect: 
		* (row, col) - - (row, col+1)
		* (row, col) - - (row+1, col)
		* (row, col) - - (row+1, col+1)
		* (row, col+1) - - (row+1, col)
		* Connect boundary:
		* if (j < ncols - 1): contains the process of the bottom row(row = nrows - 1)
		* if (i < nrows-1): contains the process of the rightmost column(col = ncols - 1)
		*/

		for (std::size_t i = 0; i < nrows; ++i) {
			for (std::size_t j = 0; j < ncols; ++j) {
				if (j < ncols - 1)set_virtual_spring(get_particle(i, j), get_particle(i, j + 1));
				if (i < nrows - 1)set_virtual_spring(get_particle(i, j), get_particle(i + 1, j));
				if ((i < nrows - 1) && (j < ncols - 1)) {
					set_virtual_spring(get_particle(i, j), get_particle(i + 1, j + 1));
					set_virtual_spring(get_particle(i, j + 1), get_particle(i + 1, j));
				}
			}
		}


		/*
		* Set the virtual spring connection between secondary neighbours
		* Does the "flexion spring" have a significant impact on cloth?
		* possible numbers of neighbours of one particle: 
		* 3(min) - - a 2*2 grid
		* 4,5,6,8,10,11,13 ... (any more values?)
		* 16(max) - - 8 direct springs and 8 second springs
		*/
		for (std::size_t i = 0; i < nrows; ++i) {
			for (std::size_t j = 0; j < ncols; ++j) {
				if (j < ncols - 2)set_virtual_spring(get_particle(i, j), get_particle(i, j + 2));
				if (i < nrows - 2)set_virtual_spring(get_particle(i, j), get_particle(i + 2, j));
				if ((i < nrows - 2) && (j < ncols - 2)) {
					set_virtual_spring(get_particle(i, j), get_particle(i + 2, j + 2));
					set_virtual_spring(get_particle(i, j + 2), get_particle(i + 2, j));
				}
			}
		}

	}

}