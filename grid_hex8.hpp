#ifndef GRID_HEX8
#define GRID_HEX8
#include <vector>

struct GridHex8Struct
{

    // point data
    int num_point = 0;
    std::vector<int> point_id_vec;
    std::vector<double> point_pos_x_vec;
    std::vector<double> point_pos_y_vec;
    std::vector<double> point_pos_z_vec;

    // element data
    int num_element = 0;
    std::vector<int> element_id_vec;
    std::vector<int> element_p00_id_vec;
    std::vector<int> element_p01_id_vec;
    std::vector<int> element_p02_id_vec;
    std::vector<int> element_p03_id_vec;
    std::vector<int> element_p04_id_vec;
    std::vector<int> element_p05_id_vec;
    std::vector<int> element_p06_id_vec;
    std::vector<int> element_p07_id_vec;

};

struct BoundaryConfigHex8Struct
{

    std::string boundary_type_str;
    std::vector<double> boundary_parameter_vec;

};

struct BoundaryHex8Struct
{

    // flux boundary condition data
    int num_element_flux = 0;
    std::vector<int> element_flux_id_vec;
    std::vector<int> element_flux_pa_loc_vec;
    std::vector<int> element_flux_pb_loc_vec;
    std::vector<int> element_flux_pc_loc_vec;
    std::vector<int> element_flux_pd_loc_vec;
    std::vector<int> element_flux_config_id_vec;
    
    // value boundary condition data
    int num_element_value = 0;
    std::vector<int> element_value_id_vec;
    std::vector<int> element_value_pa_loc_vec;
    std::vector<int> element_value_pb_loc_vec;
    std::vector<int> element_value_pc_loc_vec;
    std::vector<int> element_value_pd_loc_vec;
    std::vector<int> element_value_config_id_vec;

    // boundary condition data
    std::vector<BoundaryConfigHex8Struct> boundary_config_vec;

};

#endif
