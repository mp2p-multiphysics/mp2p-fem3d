#ifndef DOMAIN_2D
#define DOMAIN_2D
#include <fstream>
#include <set>
#include <sstream>
#include <unordered_map>
#include "container_typedef.hpp"
#include "domain_3d.hpp"

namespace FEM3D
{

class Domain2D
{
    /*

    Representation of a 2D domain.

    Variables
    =========
    domain_in : Domain3D
        3D domain where this 2D domain can be found.
    file_in_element_str_in : string
        Path to CSV file with data for domain elements.

    Notes
    ====
    The CSV file with element data must have the following columns:
        global element ID
        global point ID of local point 0
        global point ID of local point 1
        ...

    Point 0 and 1 refer to the left and right points of a line2 element.
    Point 0, 1, and 2 refer to the left, middle, and right points of a line3 element.

    */

    public:

    // 3D domain where this 2D domain is found
    Domain3D* domain_ptr;

    // type of element
    // 0 - tri3; 1 - quad4; 2 - tri6; 3 - quad8
    int type_element = 0;
    int num_neighbor = 0;

    // point data
    int num_point = 0;
    VectorInt point_pdid_to_pgid_vec;
    MapIntInt point_pgid_to_pdid_map;
    VectorDouble point_position_x_vec;
    VectorDouble point_position_y_vec;
    VectorDouble point_position_z_vec;

    // element data
    int num_element = 0;
    VectorInt element_edid_to_egid_vec;
    MapIntInt element_egid_to_edid_map;
    VectorInt2D element_edid_plid_to_pgid_vec;  // [edid][plid] -> pgid

    // file names
    std::string file_in_element_str;

    // default constructor
    Domain2D() {}

    // constructor
    Domain2D(Domain3D &domain_in, std::string file_in_element_str_in)
    {

        // store variables
        domain_ptr = &domain_in;
        file_in_element_str = file_in_element_str_in;

        // get element type
        type_element = domain_ptr->type_element;

        // get number of neighboring points per element
        switch (type_element)
        {
            case 0: num_neighbor = 3; break;  // tri3
            case 1: num_neighbor = 4; break;  // quad4
            case 2: num_neighbor = 6; break;  // tri6
            case 3: num_neighbor = 8; break;  // quad8
        }

        // get points and elements
        read_domain_element();
        extract_domain_point();

    }
    
    private:

    // functions
    void read_domain_element();
    void extract_domain_point();

};

void Domain2D::read_domain_element()
{

    // read file with elements
    std::ifstream file_in_element_stream(file_in_element_str);  

    // initialize for iteration
    bool is_element_header = true;  // true while reading header
    std::string line_element_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_element_stream, line_element_str))
    {

        // skip header
        if (is_element_header)
        {
            is_element_header = false; // not reading header
            continue;
        }

        // count number of elements
        num_element++;

        // convert line string into stringstream
        std::stringstream line_element_stream(line_element_str);

        // initialize for iteration
        int value_element_num = 0;  // counts position of value
        std::string value_element_str;  // stores values in lines
        VectorInt pgid_vec_sub;  // vector of global point IDs

        // iterate through each value
        while (std::getline(line_element_stream, value_element_str, ','))
        {

            // store values in appropriate vector
            if (value_element_num == 0)
            {
                element_edid_to_egid_vec.push_back(std::stoi(value_element_str));
            }
            else if (value_element_num < num_neighbor + 1)
            {
                pgid_vec_sub.push_back(std::stoi(value_element_str));
            }

            // store vector of global point IDs
            if (value_element_num == num_neighbor)
            {
                element_edid_plid_to_pgid_vec.push_back(pgid_vec_sub);
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    // generate map of global to domain ID for elements
    for (int edid = 0; edid < num_element; edid++)
    {
        int egid = element_edid_to_egid_vec[edid];
        element_egid_to_edid_map[egid] = edid;
    }

}

void Domain2D::extract_domain_point()
{

    // initialize set of pgid
    std::set<int> pgid_set;

    // iterate through each element and get pgid
    for (auto pgid_sub : element_edid_plid_to_pgid_vec)
    {
        pgid_set.insert(pgid_sub.begin(), pgid_sub.end());
    }

    // convert pgid set to vector
    point_pdid_to_pgid_vec = VectorInt(pgid_set.begin(), pgid_set.end());

    // iterate through each pgid and get properties
    for (int pdid = 0; pdid < point_pdid_to_pgid_vec.size(); pdid++)
    {
        
        // get pgid
        int pgid = point_pdid_to_pgid_vec[pdid];

        // get properties from 3d domain
        int pdid_3d = domain_ptr->point_pgid_to_pdid_map[pgid];
        double position_x = domain_ptr->point_position_x_vec[pdid_3d];
        double position_y = domain_ptr->point_position_y_vec[pdid_3d];
        double position_z = domain_ptr->point_position_z_vec[pdid_3d];

        // append to vectors and maps
        point_pgid_to_pdid_map[pgid] = pdid;
        point_position_x_vec.push_back(position_x);
        point_position_y_vec.push_back(position_y);
        point_position_z_vec.push_back(position_z);

    }

    // get number of points
    num_point = point_pdid_to_pgid_vec.size();

}

}

#endif
