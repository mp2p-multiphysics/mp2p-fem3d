#ifndef INPUT_GRID_HEX8_CSV
#define INPUT_GRID_HEX8_CSV
#include <sstream>
#include <fstream>
#include <vector>
#include "grid_hex8.hpp"

GridHex8Struct input_grid_hex8_csv(std::string file_in_point_str, std::string file_in_element_str)
{

    // read file with points
    std::ifstream file_in_point_stream(file_in_point_str);

    // initialize struct with grid data
    GridHex8Struct gh8s;

    // initialize for iteration
    bool is_point_header = true;  // true while reading header
    std::string line_point_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_point_stream, line_point_str))
    {

        // skip header
        if (is_point_header)
        {
            is_point_header = false; // not reading header
            continue;
        }

        // count number of particles
        gh8s.num_point++;

        // convert line string into stringstream
        std::stringstream line_point_stream(line_point_str);

        // initialize for iteration
        int value_point_num = 0;  // counts position of value
        std::string value_point_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_point_stream, value_point_str, ','))
        {

            // store values in appropriate vector
            switch (value_point_num)
            {
                case 0: gh8s.point_id_vec.push_back(std::stoi(value_point_str)); break;
                case 1: gh8s.point_pos_x_vec.push_back(std::stod(value_point_str)); break;
                case 2: gh8s.point_pos_y_vec.push_back(std::stod(value_point_str)); break;
                case 3: gh8s.point_pos_z_vec.push_back(std::stod(value_point_str)); break;
            }

            // increment value count
            value_point_num++;

        }

    }

    // close point file
    file_in_point_stream.close();

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

        // count number of particles
        gh8s.num_element++;

        // convert line string into stringstream
        std::stringstream line_element_stream(line_element_str);

        // initialize for iteration
        int value_element_num = 0;  // counts position of value
        std::string value_element_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_element_stream, value_element_str, ','))
        {

            // store values in appropriate vector
            switch (value_element_num)
            {
                case 0: gh8s.element_id_vec.push_back(std::stoi(value_element_str)); break;
                case 1: gh8s.element_p00_id_vec.push_back(std::stod(value_element_str)); break;
                case 2: gh8s.element_p01_id_vec.push_back(std::stod(value_element_str)); break;
                case 3: gh8s.element_p02_id_vec.push_back(std::stod(value_element_str)); break;
                case 4: gh8s.element_p03_id_vec.push_back(std::stod(value_element_str)); break;
                case 5: gh8s.element_p04_id_vec.push_back(std::stod(value_element_str)); break;
                case 6: gh8s.element_p05_id_vec.push_back(std::stod(value_element_str)); break;
                case 7: gh8s.element_p06_id_vec.push_back(std::stod(value_element_str)); break;
                case 8: gh8s.element_p07_id_vec.push_back(std::stod(value_element_str)); break;
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    return gh8s;

}

#endif
