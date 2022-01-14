/********************************************************************
* 
* The idea of CSF is from this paper:
* http://www.mdpi.com/2072-4292/8/6/501/htm
* The basic architectures and ideas are inspired by this article (its related code is open source).
* 
*********************************************************************/


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

// -- tin Dependencies
#include <queue>
#include <CGAL/Plane_3.h>
#include <CGAL/convex_hull_2.h>

#define M_PI 3.14159265358979323846

void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams) {

    typedef CGAL::Projection_traits_xy_3<Kernel> Gt;
    typedef CGAL::Delaunay_triangulation_2<Gt> DT;

    //build initial grid
    std::vector<const Point*> grid;
    std::vector<Point> init_tin;
    std::vector<int> class_labels;
    std::vector<Point> tin;

    double resolution = jparams["resolution"];
    double distance = jparams["distance"];
    double initial_angle = jparams["angle"];
    std::string output_las = jparams["output_las"];

    // convert angle to radian
    double angle = (initial_angle * M_PI) / 180;
    // for counting the number of ground point and non ground point
    int ground_point = 0;
    int non_ground_point = 0;
    // calculate the bounding box to create the grid
    double x_min = pointcloud[0][0];
    double y_min = pointcloud[0][1];
    double x_max = pointcloud[0][0];
    double y_max = pointcloud[0][1];
    for (int i = 1; i < pointcloud.size(); i++) {
        if (x_min > pointcloud[i][0]) x_min = pointcloud[i][0];
        if (y_min > pointcloud[i][1]) y_min = pointcloud[i][1];
        if (x_max < pointcloud[i][0]) x_max = pointcloud[i][0];
        if (y_max < pointcloud[i][1]) y_max = pointcloud[i][1];
    }
    // calculate the row and colum of the grid
    int ncols = ceil((x_max - x_min) / resolution);
    int nrows = ceil((y_max - y_min) / resolution);


    //initialize pointer to store
    for (int i = 0; i < nrows * ncols; ++i) { grid.push_back(nullptr); }

    for (int i = 0; i < pointcloud.size(); ++i) {
        int row = floor((pointcloud[i].x() - x_min) / resolution);
        int col = floor((pointcloud[i].y() - y_min) / resolution);
        int index = col + row * nrows;
        if (grid[index] == 0) grid[index] = &pointcloud[i];
        else if (grid[index] != 0) {
            if (pointcloud[i].z() < grid[index]->z()) grid[index] = &pointcloud[i];
        }
    }
    // insert initial ground point and construct DT
    DT dt;
    for (int i = 0; i < grid.size(); ++i) {
        if (grid[i] == nullptr) continue;
        else dt.insert(*grid[i]);
    }
    //calcuate the convex hull, and add the convex hull point into the DT
    std::vector<Point> convexhull_point;
    CGAL::convex_hull_2(pointcloud.begin(), pointcloud.end(), std::back_inserter(convexhull_point), Gt());
    //modify the z value
    for (int i = 0; i < convexhull_point.size(); i++) {
        int row = floor((convexhull_point[i].x() - x_min) / resolution);
        int col = floor((convexhull_point[i].y() - y_min) / resolution);
        double temp_z = grid[col + row * nrows]->z();
        Point temp = Point(convexhull_point[i][0], convexhull_point[i][1], temp_z);
        dt.insert(temp);
    }
    //  iterate every point in the point cloud data
    for (int i = 0; i < pointcloud.size(); ++i) {
        DT::Face_handle triangle = dt.locate(pointcloud[i]);
        DT::Vertex_handle v0 = triangle->vertex(0);
        DT::Vertex_handle v1 = triangle->vertex(1);
        DT::Vertex_handle v2 = triangle->vertex(2);
        double d0 = CGAL::squared_distance(v0->point(), pointcloud[i]);
        double d1 = CGAL::squared_distance(v1->point(), pointcloud[i]);
        double d2 = CGAL::squared_distance(v2->point(), pointcloud[i]);

        Kernel::Plane_3 plane_1 = Kernel::Plane_3(v0->point(), v1->point(), v2->point());
        double h = CGAL::squared_distance(pointcloud[i], plane_1);

        double temp_angle[3];
        temp_angle[0] = asin(h / d0);
        temp_angle[1] = asin(h / d1);
        temp_angle[2] = asin(h / d2);
        double max_angle = temp_angle[0];

        for (int j = 0; j < 3; j++) {
            if (temp_angle[j] > max_angle) { max_angle = temp_angle[j]; }
        }
        if (h <= distance && max_angle <= angle) {
            tin.emplace_back(pointcloud[i]);
            class_labels.emplace_back(1);
            dt.insert(pointcloud[i]);
            ground_point++;
        }
        else {
            tin.emplace_back(pointcloud[i]);
            class_labels.emplace_back(2);
            non_ground_point++;
        }
    }
    write_lasfile(jparams["output_las"], tin, class_labels);
    std::cout << "the number of ground point: " << ground_point << std::endl;
    std::cout << "there number of non ground point: " << non_ground_point << std::endl;
    std::cout << "=== TIN refinement groundfilter has finished===\n" << tin.size() << std::endl;
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
    double resolution(jparams["resolution"]);
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
    double epsilon_ground(jparams["epsilon_ground"]);
    csf::filter_classification(inverse_pointcloud, c1, epsilon_ground, class_labels);
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

