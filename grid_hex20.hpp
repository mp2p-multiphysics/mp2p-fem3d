#ifndef GRID_HEX20
#define GRID_HEX20
#include <vector>

struct GridHex20Struct
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
    std::vector<int> element_p08_id_vec;
    std::vector<int> element_p09_id_vec;
    std::vector<int> element_p10_id_vec;
    std::vector<int> element_p11_id_vec;
    std::vector<int> element_p12_id_vec;
    std::vector<int> element_p13_id_vec;
    std::vector<int> element_p14_id_vec;
    std::vector<int> element_p15_id_vec;
    std::vector<int> element_p16_id_vec;
    std::vector<int> element_p17_id_vec;
    std::vector<int> element_p18_id_vec;
    std::vector<int> element_p19_id_vec;

};

struct BoundaryConfigHex20Struct
{

    std::string boundary_type_str;
    std::vector<double> boundary_parameter_vec;

};

struct BoundaryHex20Struct
{

    // flux boundary condition data
    int num_element_flux = 0;
    std::vector<int> element_flux_id_vec;
    std::vector<int> element_flux_pa_loc_vec;
    std::vector<int> element_flux_pb_loc_vec;
    std::vector<int> element_flux_pc_loc_vec;
    std::vector<int> element_flux_pd_loc_vec;
    std::vector<int> element_flux_pe_loc_vec;
    std::vector<int> element_flux_pf_loc_vec;
    std::vector<int> element_flux_pg_loc_vec;
    std::vector<int> element_flux_ph_loc_vec;
    std::vector<int> element_flux_config_id_vec;
    
    // value boundary condition data
    int num_element_value = 0;
    std::vector<int> element_value_id_vec;
    std::vector<int> element_value_pa_loc_vec;
    std::vector<int> element_value_pb_loc_vec;
    std::vector<int> element_value_pc_loc_vec;
    std::vector<int> element_value_pd_loc_vec;
    std::vector<int> element_value_pe_loc_vec;
    std::vector<int> element_value_pf_loc_vec;
    std::vector<int> element_value_pg_loc_vec;
    std::vector<int> element_value_ph_loc_vec;
    std::vector<int> element_value_config_id_vec;

    // boundary condition data
    std::vector<BoundaryConfigHex20Struct> boundary_config_vec;

};

#endif
