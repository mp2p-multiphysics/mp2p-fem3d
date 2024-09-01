#ifndef INPUT_GRID_HEX20_CSV
#define INPUT_GRID_HEX20_CSV
#include <sstream>
#include <fstream>
#include <vector>
#include "grid_hex20.hpp"

GridHex20Struct input_grid_hex20_csv(std::string file_in_point_str, std::string file_in_element_str)
{

    // read file with points
    std::ifstream file_in_point_stream(file_in_point_str);

    // initialize struct with grid data
    GridHex20Struct gh20s;

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
        gh20s.num_point++;

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
                case 0: gh20s.point_id_vec.push_back(std::stoi(value_point_str)); break;
                case 1: gh20s.point_pos_x_vec.push_back(std::stod(value_point_str)); break;
                case 2: gh20s.point_pos_y_vec.push_back(std::stod(value_point_str)); break;
                case 3: gh20s.point_pos_z_vec.push_back(std::stod(value_point_str)); break;
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
        gh20s.num_element++;

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
                case 0: gh20s.element_id_vec.push_back(std::stoi(value_element_str)); break;
                case 1: gh20s.element_p00_id_vec.push_back(std::stod(value_element_str)); break;
                case 2: gh20s.element_p01_id_vec.push_back(std::stod(value_element_str)); break;
                case 3: gh20s.element_p02_id_vec.push_back(std::stod(value_element_str)); break;
                case 4: gh20s.element_p03_id_vec.push_back(std::stod(value_element_str)); break;
                case 5: gh20s.element_p04_id_vec.push_back(std::stod(value_element_str)); break;
                case 6: gh20s.element_p05_id_vec.push_back(std::stod(value_element_str)); break;
                case 7: gh20s.element_p06_id_vec.push_back(std::stod(value_element_str)); break;
                case 8: gh20s.element_p07_id_vec.push_back(std::stod(value_element_str)); break;
                case 9: gh20s.element_p08_id_vec.push_back(std::stod(value_element_str)); break;
                case 10: gh20s.element_p09_id_vec.push_back(std::stod(value_element_str)); break;
                case 11: gh20s.element_p10_id_vec.push_back(std::stod(value_element_str)); break;
                case 12: gh20s.element_p11_id_vec.push_back(std::stod(value_element_str)); break;
                case 13: gh20s.element_p12_id_vec.push_back(std::stod(value_element_str)); break;
                case 14: gh20s.element_p13_id_vec.push_back(std::stod(value_element_str)); break;
                case 15: gh20s.element_p14_id_vec.push_back(std::stod(value_element_str)); break;
                case 16: gh20s.element_p15_id_vec.push_back(std::stod(value_element_str)); break;
                case 17: gh20s.element_p16_id_vec.push_back(std::stod(value_element_str)); break;
                case 18: gh20s.element_p17_id_vec.push_back(std::stod(value_element_str)); break;
                case 19: gh20s.element_p18_id_vec.push_back(std::stod(value_element_str)); break;
                case 20: gh20s.element_p19_id_vec.push_back(std::stod(value_element_str)); break;                
            }

            // increment value count
            value_element_num++;

        }

    }

    // close element file
    file_in_element_stream.close();

    return gh20s;

}

#endif
