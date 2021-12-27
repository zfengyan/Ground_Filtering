/*
  GEO1015.2021
  hw03 
  --
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
*/

#include "GroundFilter.h"

// -- LAS reading and writing
#include <lasreader.hpp>
#include <laswriter.hpp>

// -- CGAL delaunay triangulation
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

// -- CGAL kd-tree
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Search_traits_2.h>

// -- csf Dependencies
#include "Cloth.h"
#include <cmath>


void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
    !!! TO BE COMPLETED !!!
      
    Function that performs ground filtering using TIN refinement and writes the result to a new LAS file.

    !!! You are free to subdivide the functionality of this function into several functions !!!
      
    Inputs:
      pointcloud: input point cloud (an Nx3 numpy array),
      jparams: a dictionary jparams with all the parameters that are to be used in this function:
        - resolution:    resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
        - distance:      distance threshold used in the ground filtering algorithm,
        - angle:         angle threshold used in the ground filtering algorithm in degrees,
        - output_las:    path to output .las file that contains your ground classification,
  */
  typedef CGAL::Projection_traits_xy_3<Kernel>  Gt;
  typedef CGAL::Delaunay_triangulation_2<Gt> DT;

  // double resolution = j["resolution"];
  // double distance = j["distance"];
  // double angle = j["angle"];
  // std::string output_las = jparams["output_las"];

  //-- TIP CGAL triangulation -> https://doc.cgal.org/latest/Triangulation_2/index.html
  //-- Insert points in a triangulation: [https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a1025cd7e7226ccb44d82f0fb1d63ad4e]
  // DT dt;
  // dt.insert(Point(0,0,0));
  // dt.insert(Point(10,0,0));
  // dt.insert(Point(0,10,0));
  //-- Find triangle that intersects a given point: [https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a940567120751e7864c7b345eaf756642]
  // DT::Face_handle triangle = dt.locate(Point(3,3,0));
  //-- get the 3 vertices of the triangle: 
  // DT::Vertex_handle v0 = triangle->vertex(0);
  // DT::Vertex_handle v1 = triangle->vertex(1);
  // DT::Vertex_handle v2 = triangle->vertex(2);
  // get the coordinates of the three vertices:
  // std::cout << "v0 has the coordinates ( " << v0->point().x() << "  " << v0->point().y() << " " << v0->point().z() << " )" << std::endl;
  // std::cout << "v1 has the coordinates ( " << v1->point().x() << "  " << v1->point().y() << " " << v1->point().z() << " )" << std::endl;
  // std::cout << "v2 has the coordinates ( " << v2->point().x() << "  " << v2->point().y() << " " << v2->point().z() << " )" << std::endl;

  //-- TIP CGAL compute squared distance between two points: [https://doc.cgal.org/latest/Kernel_23/group__squared__distance__grp.html#ga1ff73525660a052564d33fbdd61a4f71]
  // std::cout << "the squared distance between v0 and v1 is: " << CGAL::squared_distance(v0->point(), v1->point()) << std::endl;
  
  //-- TIP
  //-- write the results to a new LAS file
  // std::vector<int> class_labels;
  // write_lasfile(jparams["output_las"], pointcloud, class_labels);
}


// related functions of CSF algorithm
namespace csf {
    
    /*
    * @brief: compute the bounding box points("low left" corner and "up right" corner in 3D space)
    * @param: pcloud: input point cloud (an N x 3 array)
    *         bmin: "low-left" corner bmax: "up-right" corner
    * @return: none
    */
    void bounding_box(const std::vector<Point>& pcloud, MyPoint& bmin, MyPoint& bmax) {

        if (!pcloud.empty()) {
            // bmin = bmax = pcloud[0];
            bmin.x = bmax.x = pcloud[0][0]; // pcloud[0][0]: x of point pcloud[0]
            bmin.y = bmax.y = pcloud[0][1]; // y
            bmin.z = bmax.z = pcloud[0][2]; // z

            for (std::size_t i = 1; i < pcloud.size(); ++i) {
                const Point& tmp = pcloud[i]; // reference point at the beginning
                if (pcloud[i][0] < bmin.x)bmin.x = pcloud[i][0];
                else if(pcloud[i][0] > bmax.x)bmax.x = pcloud[i][0];

                if (pcloud[i][1] < bmin.y)bmin.y = pcloud[i][1];
                else if (pcloud[i][1] > bmax.y)bmax.y = pcloud[i][1];

                if (pcloud[i][2] < bmin.z)bmin.z = pcloud[i][2];
                else if(pcloud[i][2] > bmax.z)bmax.z = pcloud[i][2];
            }
        }

    }


    void write_lasfile_tmp(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels) {
        /*
        Function to write a new LAS file with point labels (for the LAS classification field)

        Inputs:
          filename:   the filename to write the LAS file to
          pointcloud: input point cloud (a vector of Points),
          Labels:     Contains point labels. Should be a vector of ints of the same size as pointcloud (ie. one label for each point in the same order as pointcloud). Uses LAS classification codes, ie 2 = ground. 1 = unclassified.
        */
        LASwriteOpener laswriteopener;
        laswriteopener.set_file_name(filename.c_str());

        LASheader lasheader;
        lasheader.x_scale_factor = 0.01;
        lasheader.y_scale_factor = 0.01;
        lasheader.z_scale_factor = 0.01;
        lasheader.x_offset = 0.0;
        lasheader.y_offset = 0.0;
        lasheader.z_offset = 0.0;
        lasheader.point_data_format = 0;
        lasheader.point_data_record_length = 20;

        LASpoint laspoint;
        laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

        LASwriter* laswriter = laswriteopener.open(&lasheader);
        if (laswriter == 0)
        {
            std::cerr << "ERROR: could not open laswriter\n";
            exit(1);
        }

        if (pointcloud.size() != class_labels.size()) {
            std::cerr << "ERROR: points has a different size than class_labels\n";
            exit(1);
        }

        for (size_t i = 0; i < pointcloud.size(); ++i) {
            const Point& p = pointcloud[i];
            const int& label = class_labels[i];

            laspoint.set_x(p[0]);
            laspoint.set_y(p[1]);
            laspoint.set_z(p[2]);

            laspoint.set_classification(label);

            laswriter->write_point(&laspoint);
            laswriter->update_inventory(&laspoint);
        }

        laswriter->update_header(&lasheader, TRUE);
        laswriter->close();
        delete laswriter;
    }


    /*
    * output a cloth(its particles)
    */
    void write_lasfile_particles(const std::string filename, const std::vector<Particle>& particles, const std::vector<int>& class_labels) {
     
        LASwriteOpener laswriteopener;
        laswriteopener.set_file_name(filename.c_str());

        LASheader lasheader;
        lasheader.x_scale_factor = 0.01;
        lasheader.y_scale_factor = 0.01;
        lasheader.z_scale_factor = 0.01;
        lasheader.x_offset = 0.0;
        lasheader.y_offset = 0.0;
        lasheader.z_offset = 0.0;
        lasheader.point_data_format = 0;
        lasheader.point_data_record_length = 20;

        LASpoint laspoint;
        laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

        LASwriter* laswriter = laswriteopener.open(&lasheader);
        if (laswriter == 0)
        {
            std::cerr << "ERROR: could not open laswriter\n";
            exit(1);
        }

        if (particles.size() != class_labels.size()) {
            std::cerr << "ERROR: points has a different size than class_labels\n";
            exit(1);
        }

        for (size_t i = 0; i < particles.size(); ++i) {
            const Particle& p = particles[i];
            const int& label = class_labels[i];

            laspoint.set_x(p.cur_pos.v[0]);
            laspoint.set_y(p.cur_pos.v[1]);
            laspoint.set_z(p.cur_pos.v[2]);
            //laspoint.set_z(p.Intersect_Height_Value);

            laspoint.set_classification(label);

            laswriter->write_point(&laspoint);
            laswriter->update_inventory(&laspoint);
        }

        laswriter->update_header(&lasheader, TRUE);
        laswriter->close();
        delete laswriter;
    }


    /*
    * @brief: 
    * Projecting all the LiDAR points and grid particles to a horizontal plane 
    * and finding the CP(corresponding point) for each grid particle in this plane. 
    * then recording the Intersection Height Value.
    * @param: pointcloud, cloth
    */
    void find_intersection_height(const std::vector<Point>& pointcloud, Cloth& cloth) {

        typedef CGAL::Simple_cartesian<double> K;
        typedef K::Point_2 Point_2;

        typedef boost::tuple<Point_2, std::size_t> Point_and_int;
        typedef CGAL::Search_traits_2<K> Traits_base; //typedef CGAL::Search_traits_3<Kernel> Traits_base;     
        typedef CGAL::Search_traits_adapter<Point_and_int,
            CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
            Traits_base> Traits;
        typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
        typedef Neighbor_search::Tree Tree;


        /*typedef CGAL::Search_traits_3<Kernel> TreeTraits;
        typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
        typedef Neighbor_search::Tree Tree;*/

        // project the pointcloud to a plane: Point_2 type
        std::vector<Point_2> pointcloud_proj;
        pointcloud_proj.reserve(pointcloud.size());
        for (auto point : pointcloud) {
            pointcloud_proj.emplace_back(Point_2(point[0], point[1]));
        }
       

        // set the indices vector
        std::vector<std::size_t>indices;
        indices.reserve(pointcloud.size());
        for (std::size_t i = 0; i < pointcloud.size(); ++i) {
            indices.emplace_back(i);
        }


        //Tree tree(pointcloud.begin(), pointcloud.end());
        Tree tree(boost::make_zip_iterator(boost::make_tuple(pointcloud_proj.begin(), indices.begin())),
                  boost::make_zip_iterator(boost::make_tuple(pointcloud_proj.end(), indices.end())));

        const unsigned int N = 1;
        
        for (std::size_t i = 0; i < cloth.particles.size(); ++i) {

            // query_point: projection of particle
            Point_2 query_point(cloth.particles[i].cur_pos.v[0], cloth.particles[i].cur_pos.v[1]); //x, y
            Neighbor_search search_result(tree, query_point, N);
            auto res(search_result.begin()); // N=1, only one neighbour
            std::size_t indice(boost::get<1>(res->first)); // indice of the result point
            cloth.particles[i].Intersect_Height_Value = pointcloud[indice][2]; // set the intersection height value
        }

        //for (auto res : search_result) {
        //    /*Point neighbour_point(res.first);
        //    double distance(res.second);
        //    std::cout << "neighbouring distance: " << res.first<< '\n';*/
        //}

        /*for (auto it = search_result.begin(); it != search_result.end(); ++it) {
            std::cout << " d(q, nearest neighbor)=  "
                << boost::get<0>(it->first) << " " << boost::get<1>(it->first) << std::endl;
        }*/
    }


    /*
    * @brief: compute the distance between two pointclouds
    * compare the distance to the threshold
    * store the classification results in the vector
    * '2' for ground points
    * '1' for other points
    */
    void filter_classification(
        const std::vector<Point>& inverse_pointcloud,
        Cloth& cloth,
        const double& threshold,
        std::vector<int>& class_labels) {

        class_labels.reserve(inverse_pointcloud.size());

        for (std::size_t i = 0; i < inverse_pointcloud.size(); ++i) {

            // find the x and y shift from the initial_position of the cloth
            double delta_x(inverse_pointcloud[i][0] - cloth.initial_position.v[0]);
            double delta_y(inverse_pointcloud[i][1] - cloth.initial_position.v[1]);

            // find the approximate position in the cloth frid
            int row(int(delta_x / cloth.row_step));
            int col(int(delta_y / cloth.col_step));

            // neighboring point in cloth
            int row_right(row), col_right(col + 1); // right neighbor
            int row_diagonal(row + 1), col_diagonal(col + 1); // diagonal neighbor
            int row_bottom(row + 1), col_bottom(col); //bottom neighbor

            // calculate the "residuals"
            double residual_x((delta_x - row * cloth.row_step) / cloth.row_step);
            double residual_y((delta_y - col * cloth.col_step) / cloth.col_step);

            double distance(
                cloth.get_particle(row, col)->cur_pos.v[2] * (1 - residual_x) * (1 - residual_y) +
                cloth.get_particle(row_right, col_right)->cur_pos.v[2] * residual_x * (1 - residual_y) +
                cloth.get_particle(row_diagonal, col_diagonal)->cur_pos.v[2] * residual_x * residual_y +
                cloth.get_particle(row_bottom, col_bottom)->cur_pos.v[2] * (1 - residual_x) * residual_y
            );
            double variance(std::fabs(distance - inverse_pointcloud[i][2]));

            // classification
            if (variance < threshold) { class_labels.emplace_back(2); }
            else { class_labels.emplace_back(1); }
        }
    }

}


void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams) {
    /*
    !!! TO BE COMPLETED !!!
    
    Function that performs ground filtering using CSF and writes the result to a new LAS file.

    !!! You are free to subdivide the functionality of this function into several functions !!!
    
    Inputs:
        pointcloud: input point cloud (an Nx3 numpy array),
        jparams: a dictionary with all the parameters that are to be used in this function:
        - resolution:     resolution of the cloth grid,
        - epsilon_zmax:   tolerance to stop the iterations,
        - epsilon_ground: threshold used to classify ground points,
        - output_las:     path to output .las file that contains your ground classification
    */

    /*typedef CGAL::Search_traits_3<Kernel> TreeTraits;
    typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef Neighbor_search::Tree Tree;*/

    // double resolution = j["resolution"];
    // double epsilon_zmax = j["epsilon_zmax"];
    // double epsilon_ground = j["epsilon_ground"];
    // std::string output_las = jparams["output_las"];

    // //-- print the first 5 points in the pointcloud, which are CGAL Point_3
    // //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html


    //-- TIP
    //-- construct and query kd-tree:
    // https://doc.cgal.org/latest/Spatial_searching/index.html#title5
    //

    /*
    * Inverse the original pointcloud and output
    */
    std::vector<Point> inverse_pointcloud;
    inverse_pointcloud.reserve(pointcloud.size());
    std::size_t i = 0;
    for (auto p : pointcloud) {
        inverse_pointcloud.emplace_back(Point(p[0], p[1], -p[2]));
        ++i;
        if (i == pointcloud.size())
            break;
    }

    csf::MyPoint pmin, pmax;
    bounding_box(inverse_pointcloud, pmin, pmax);
    double resolution(0.5);
    std::size_t NROWS((std::size_t)(ceil((pmax.x - pmin.x) / resolution) + 1));
    std::size_t NCOLS((std::size_t)(ceil((pmax.y - pmin.y) / resolution) + 1));
    csf::Cloth c1(NROWS, NCOLS, 1, resolution, resolution, 0.01, csf::Vector3d(pmin.x, pmin.y, pmax.z));

    find_intersection_height(inverse_pointcloud, c1);
    //add gravity
    c1.addforce_for_particles(csf::Vector3d(0, 0, -10));

    //int count(0);
    for (std::size_t i = 0; i < 500; ++i) {
        c1.update_cloth_gravity();
        c1.update_cloth_spring();
        c1.terrain_intersection_check();
        
        //if (c1.calculate_max_diff() != 0 && c1.calculate_max_diff() >= 0.066)break;
        //++count;
    }
    //std::cout << c1.calculate_max_diff();


    std::vector<int> class_labels;
    csf::filter_classification(inverse_pointcloud, c1, 0.5, class_labels);
    write_lasfile(jparams["output_las"], pointcloud, class_labels);

    // output inverse_cloud
    //std::vector<int> class_labels(inverse_pointcloud.size()); // Initialized with 0
    //csf::write_lasfile_tmp(jparams["output_las"], inverse_pointcloud, class_labels);

    // output the simulated cloth
    //std::vector<int> class_labels(c1.particles.size()); // Initialized with 0
    //csf::write_lasfile_particles(jparams["output_las"], c1.particles, class_labels);

 

    //-- TIP
    //-- write the results to a new LAS file
    //std::vector<int> class_labels(pointcloud.size()); // Initialized with 0
    //write_lasfile(jparams["output_las"], pointcloud, class_labels);


}



std::vector<Point> read_lasfile(const json& jparams) {
  /*
  Function to read points from a LAS file

  Inputs:
    jparams["filename"]:   the filename to read the LAS file to

  Returns:
    a std::vector<Point> with the points from the LAS file
  */
  std::string filename = jparams["filename"];
	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(filename.c_str());
	LASreader* lasreader = lasreadopener.open();
	
	if (!lasreader){
		std::cerr << "cannot read las file: " << filename << "\n";
		exit(1);
	}

  //-- store each point in a CGAL Point_3 object
  //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
	std::vector<Point> points;
	while (lasreader->read_point()) {
		points.emplace_back( 
			Point(
				lasreader->point.get_x(),
				lasreader->point.get_y(),
				lasreader->point.get_z()
			)
		);
	}
	lasreader->close();
	delete lasreader;

	return points;
}



void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels) {
  /*
  Function to write a new LAS file with point labels (for the LAS classification field)

  Inputs:
    filename:   the filename to write the LAS file to
    pointcloud: input point cloud (a vector of Points),
    Labels:     Contains point labels. Should be a vector of ints of the same size as pointcloud (ie. one label for each point in the same order as pointcloud). Uses LAS classification codes, ie 2 = ground. 1 = unclassified.
  */
  LASwriteOpener laswriteopener;
  laswriteopener.set_file_name(filename.c_str());

  LASheader lasheader;
  lasheader.x_scale_factor = 0.01;
  lasheader.y_scale_factor = 0.01;
  lasheader.z_scale_factor = 0.01;
  lasheader.x_offset = 0.0;
  lasheader.y_offset = 0.0;
  lasheader.z_offset = 0.0;
  lasheader.point_data_format = 0;
  lasheader.point_data_record_length = 20;

  LASpoint laspoint;
  laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

  LASwriter* laswriter = laswriteopener.open(&lasheader);
  if (laswriter == 0)
  {
    std::cerr << "ERROR: could not open laswriter\n";
    exit(1);
  }

	if (pointcloud.size()!=class_labels.size()) {
		std::cerr << "ERROR: points has a different size than class_labels\n";
		exit(1);
	}

  for (size_t i=0; i<pointcloud.size(); ++i) {
		const Point& p = pointcloud[i];
		const int& label = class_labels[i];

        laspoint.set_x(p[0]);
        laspoint.set_y(p[1]);
        laspoint.set_z(p[2]);

        /*
        laspoint.set_x(p.x);
        laspoint.set_y(p.y);
        laspoint.set_z(p.z);*/

		laspoint.set_classification(label);

        laswriter->write_point(&laspoint);
        laswriter->update_inventory(&laspoint);    
  } 

    laswriter->update_header(&lasheader, TRUE);
    laswriter->close();
    delete laswriter;
}

