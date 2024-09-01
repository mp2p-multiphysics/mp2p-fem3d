#ifndef MODULE_HEAT_STEADY_HEX8
#define MODULE_HEAT_STEADY_HEX8
#include <vector>
#include "Eigen/Eigen"
#include "integral_hex8.hpp"
#include "grid_hex8.hpp"
#include "scalar_hex8.hpp"

class HeatSteadyHex8Class
{

    public:

    // variables
    GridHex8Struct gh8s;
    BoundaryHex8Struct bh8s;
    IntegralHex8Class ih8c;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
        ScalarHex8Class &thermal_conductivity_sh8c, ScalarHex8Class &heat_generation_sh8c,
        int start_id
    );

    // constructor
    HeatSteadyHex8Class(GridHex8Struct &gh8s_in, BoundaryHex8Struct &bh8s_in, IntegralHex8Class &ih8c_in)
    {
        
        // variables
        gh8s = gh8s_in;
        bh8s = bh8s_in;
        ih8c = ih8c_in;

        // integrals
        ih8c.evaluate_test_functions_derivatives();
        ih8c.evaluate_integral_div_Ni_hex8_dot_div_Nj_hex8();
        ih8c.evaluate_integral_Ni_hex8();
        
    }

};

void HeatSteadyHex8Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
    ScalarHex8Class &thermal_conductivity_sh8c, ScalarHex8Class &heat_generation_sh8c,
    int start_id
)
{

    // iterate for each grid element
    for (int m = 0; m < gh8s.num_element; m++)
    {

        // get id of points around element
        int n00 = gh8s.element_p00_id_vec[m];
        int n01 = gh8s.element_p01_id_vec[m];
        int n02 = gh8s.element_p02_id_vec[m];
        int n03 = gh8s.element_p03_id_vec[m];
        int n04 = gh8s.element_p04_id_vec[m];
        int n05 = gh8s.element_p05_id_vec[m];
        int n06 = gh8s.element_p06_id_vec[m];
        int n07 = gh8s.element_p07_id_vec[m];
        int n_arr[8] = {n00, n01, n02, n03, n04, n05, n06, n07};

        // get thermal conductivity of points around element
        double thermcond00 = thermal_conductivity_sh8c.scalar_vec[n00];
        double thermcond01 = thermal_conductivity_sh8c.scalar_vec[n01];
        double thermcond02 = thermal_conductivity_sh8c.scalar_vec[n02];
        double thermcond03 = thermal_conductivity_sh8c.scalar_vec[n03];
        double thermcond04 = thermal_conductivity_sh8c.scalar_vec[n04];
        double thermcond05 = thermal_conductivity_sh8c.scalar_vec[n05];
        double thermcond06 = thermal_conductivity_sh8c.scalar_vec[n06];
        double thermcond07 = thermal_conductivity_sh8c.scalar_vec[n07];
        double thermcond_arr[8] = {thermcond00, thermcond01, thermcond02, thermcond03, thermcond04, thermcond05, thermcond06, thermcond07};

        // get heat generation of points around element
        double heatgen00 = heat_generation_sh8c.scalar_vec[n00];
        double heatgen01 = heat_generation_sh8c.scalar_vec[n01];
        double heatgen02 = heat_generation_sh8c.scalar_vec[n02];
        double heatgen03 = heat_generation_sh8c.scalar_vec[n03];
        double heatgen04 = heat_generation_sh8c.scalar_vec[n04];
        double heatgen05 = heat_generation_sh8c.scalar_vec[n05];
        double heatgen06 = heat_generation_sh8c.scalar_vec[n06];
        double heatgen07 = heat_generation_sh8c.scalar_vec[n07];
        double heatgen_arr[8] = {heatgen00, heatgen01, heatgen02, heatgen03, heatgen04, heatgen05, heatgen06, heatgen07};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable

        // calculate a_mat coefficients
        for (int i = 0; i < 8; i++){
        for (int j = 0; j < 8; j++){
            a_mat.coeffRef(n_arr[i], n_arr[j]) += thermcond_arr[i]*ih8c.integral_div_Ni_hex8_dot_div_Nj_hex8_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++)
        {
            b_vec.coeffRef(n_arr[i]) += heatgen_arr[i]*ih8c.integral_Ni_hex8_vec[m][i];
        }

    }

    // iterate for each flux boundary element
    for (int k = 0; k < bh8s.num_element_flux; k++)
    {

        // get id of element
        int m = bh8s.element_flux_id_vec[k];

        // get id of points around element
        int n00 = gh8s.element_p00_id_vec[m];
        int n01 = gh8s.element_p01_id_vec[m];
        int n02 = gh8s.element_p02_id_vec[m];
        int n03 = gh8s.element_p03_id_vec[m];
        int n04 = gh8s.element_p04_id_vec[m];
        int n05 = gh8s.element_p05_id_vec[m];
        int n06 = gh8s.element_p06_id_vec[m];
        int n07 = gh8s.element_p07_id_vec[m];
        int n_arr[8] = {n00, n01, n02, n03, n04, n05, n06, n07};

        // get x coordinates of points around element
        double x00 = gh8s.point_pos_x_vec[n00];
        double x01 = gh8s.point_pos_x_vec[n01];
        double x02 = gh8s.point_pos_x_vec[n02];
        double x03 = gh8s.point_pos_x_vec[n03];
        double x04 = gh8s.point_pos_x_vec[n04];
        double x05 = gh8s.point_pos_x_vec[n05];
        double x06 = gh8s.point_pos_x_vec[n06];
        double x07 = gh8s.point_pos_x_vec[n07];
        double x_arr[8] = {x00, x01, x02, x03, x04, x05, x06, x07};

        // get y coordinates of points around element
        double y00 = gh8s.point_pos_y_vec[n00];
        double y01 = gh8s.point_pos_y_vec[n01];
        double y02 = gh8s.point_pos_y_vec[n02];
        double y03 = gh8s.point_pos_y_vec[n03];
        double y04 = gh8s.point_pos_y_vec[n04];
        double y05 = gh8s.point_pos_y_vec[n05];
        double y06 = gh8s.point_pos_y_vec[n06];
        double y07 = gh8s.point_pos_y_vec[n07];
        double y_arr[8] = {y00, y01, y02, y03, y04, y05, y06, y07};

        // get z coordinates of points around element
        double z00 = gh8s.point_pos_z_vec[n00];
        double z01 = gh8s.point_pos_z_vec[n01];
        double z02 = gh8s.point_pos_z_vec[n02];
        double z03 = gh8s.point_pos_z_vec[n03];
        double z04 = gh8s.point_pos_z_vec[n04];
        double z05 = gh8s.point_pos_z_vec[n05];
        double z06 = gh8s.point_pos_z_vec[n06];
        double z07 = gh8s.point_pos_z_vec[n07];
        double z_arr[8] = {z00, z01, z02, z03, z04, z05, z06, z07};

        // get points where the boundary is applied
        int a = bh8s.element_flux_pa_loc_vec[k];  // 0 to 7
        int b = bh8s.element_flux_pb_loc_vec[k];  // 0 to 7
        int c = bh8s.element_flux_pc_loc_vec[k];  // 0 to 7
        int d = bh8s.element_flux_pd_loc_vec[k];  // 0 to 7

        // identify face where boundary is applied
        int local_point_int = 1000*a + 100*b + 10*c + d;
        int face_id = 0;
        switch(local_point_int)
        {
            case 0145: face_id = 0; break;  // -a
            case 2367: face_id = 1; break;  // +a
            case 0123: face_id = 2; break;  // -b
            case 4567: face_id = 3; break;  // +b
            case 0347: face_id = 4; break;  // -c
            case 1256: face_id = 5; break;  // +c
        }

        // identify boundary type
        int config_id = bh8s.element_flux_config_id_vec[k];
        BoundaryConfigHex8Struct bch8s = bh8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch8s.boundary_type_str == "neumann")
        {

            // add to b_vec
            b_vec.coeffRef(n_arr[a]) += -bch8s.boundary_parameter_vec[0] * ih8c.integral_surface_Ni_hex8(a, face_id, x_arr, y_arr, z_arr);
            b_vec.coeffRef(n_arr[b]) += -bch8s.boundary_parameter_vec[0] * ih8c.integral_surface_Ni_hex8(b, face_id, x_arr, y_arr, z_arr);
            b_vec.coeffRef(n_arr[c]) += -bch8s.boundary_parameter_vec[0] * ih8c.integral_surface_Ni_hex8(c, face_id, x_arr, y_arr, z_arr);
            b_vec.coeffRef(n_arr[d]) += -bch8s.boundary_parameter_vec[0] * ih8c.integral_surface_Ni_hex8(d, face_id, x_arr, y_arr, z_arr);

        }

    }

    // clear rows with value boundary elements
    for (int k = 0; k < bh8s.num_element_value; k++)
    {

        // get id of element
        int m = bh8s.element_value_id_vec[k];

        // get id of points around element
        int n00 = gh8s.element_p00_id_vec[m];
        int n01 = gh8s.element_p01_id_vec[m];
        int n02 = gh8s.element_p02_id_vec[m];
        int n03 = gh8s.element_p03_id_vec[m];
        int n04 = gh8s.element_p04_id_vec[m];
        int n05 = gh8s.element_p05_id_vec[m];
        int n06 = gh8s.element_p06_id_vec[m];
        int n07 = gh8s.element_p07_id_vec[m];
        int n_arr[8] = {n00, n01, n02, n03, n04, n05, n06, n07};

        // get points where the boundary is applied
        int a = bh8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = bh8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = bh8s.element_value_pc_loc_vec[k];  // 0 to 7
        int d = bh8s.element_value_pd_loc_vec[k];  // 0 to 7

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(n_arr[a]) *= 0.;
            b_vec.coeffRef(n_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(n_arr[b]) *= 0.;
            b_vec.coeffRef(n_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(n_arr[c]) *= 0.;
            b_vec.coeffRef(n_arr[c]) = 0.;
        }
        if (d != -1)
        {
            a_mat.row(n_arr[d]) *= 0.;
            b_vec.coeffRef(n_arr[d]) = 0.;
        }

    }

    // iterate for each value boundary element
    for (int k = 0; k < bh8s.num_element_value; k++)
    {

        // get id of element
        int m = bh8s.element_value_id_vec[k];

        // get id of points around element
        int n00 = gh8s.element_p00_id_vec[m];
        int n01 = gh8s.element_p01_id_vec[m];
        int n02 = gh8s.element_p02_id_vec[m];
        int n03 = gh8s.element_p03_id_vec[m];
        int n04 = gh8s.element_p04_id_vec[m];
        int n05 = gh8s.element_p05_id_vec[m];
        int n06 = gh8s.element_p06_id_vec[m];
        int n07 = gh8s.element_p07_id_vec[m];
        int n_arr[8] = {n00, n01, n02, n03, n04, n05, n06, n07};

        // get points where the boundary is applied
        int a = bh8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = bh8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = bh8s.element_value_pc_loc_vec[k];  // 0 to 7
        int d = bh8s.element_value_pd_loc_vec[k];  // 0 to 7
        
        // identify boundary type
        int config_id = bh8s.element_value_config_id_vec[k];
        BoundaryConfigHex8Struct bch8s = bh8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch8s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(n_arr[a], n_arr[a]) += 1.;
                b_vec.coeffRef(n_arr[a]) += bch8s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(n_arr[b], n_arr[b]) += 1.;
                b_vec.coeffRef(n_arr[b]) += bch8s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(n_arr[c], n_arr[c]) += 1.;
                b_vec.coeffRef(n_arr[c]) += bch8s.boundary_parameter_vec[0];
            }
            if (d != -1)
            {
                a_mat.coeffRef(n_arr[d], n_arr[d]) += 1.;
                b_vec.coeffRef(n_arr[d]) += bch8s.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
