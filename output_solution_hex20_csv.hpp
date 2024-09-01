#ifndef OUTPUT_SOLUTION_HEX20_CSV
#define OUTPUT_SOLUTION_HEX20_CSV
#include <vector>
#include <fstream>
#include "Eigen/Eigen"
#include "grid_hex20.hpp"
#include "scalar_hex20.hpp"

void output_solution_hex20_csv(ScalarHex20Class &sh20c, std::string file_out_base_str)
{

        // output file name
        std::string file_out_str = file_out_base_str + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridHex20Struct gh20s = sh20c.gh20s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,pos_z,value\n";
        for (int n = 0; n < gh20s.num_point; n++)
        {
            file_out_stream << gh20s.point_id_vec[n] << "," << gh20s.point_pos_x_vec[n] << "," << gh20s.point_pos_y_vec[n] << "," << gh20s.point_pos_z_vec[n] << "," << sh20c.scalar_vec[n] << "\n";
        }

}

void output_solution_hex20_csv(ScalarHex20Class &sh20c, std::string file_out_base_str, int ts)
{

        // output file name
        std::string file_out_str = file_out_base_str + std::to_string(ts) + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridHex20Struct gh20s = sh20c.gh20s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,pos_z,value\n";
        for (int n = 0; n < gh20s.num_point; n++)
        {
            file_out_stream << gh20s.point_id_vec[n] << "," << gh20s.point_pos_x_vec[n] << "," << gh20s.point_pos_y_vec[n] << "," << gh20s.point_pos_z_vec[n] << "," << sh20c.scalar_vec[n] << "\n";
        }

}

#endif
