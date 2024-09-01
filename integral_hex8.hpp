#ifndef INTEGRAL_HEX8
#define INTEGRAL_HEX8
#include <vector>
#include "grid_hex8.hpp"
#include "Eigen/Eigen"

class IntegralHex8Class
{

    public:
    
    // variables
    GridHex8Struct gh8s;

    // vectors with test functions and derivatives
    std::vector<std::vector<double>> jacobian_determinant_vec;
    std::vector<std::vector<std::vector<double>>> N_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_x_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_y_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_z_vec;

    // vectors with integrals
    std::vector<std::vector<double>> integral_Ni_hex8_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_hex8_Nj_hex8_vec;
    std::vector<std::vector<std::vector<double>>> integral_div_Ni_hex8_dot_div_Nj_hex8_vec;

    // functions for computing integrals
    void evaluate_test_functions_derivatives();
    void evaluate_integral_Ni_hex8();
    void evaluate_integral_Ni_hex8_Nj_hex8();
    void evaluate_integral_div_Ni_hex8_dot_div_Nj_hex8();

    // functions for computing boundary integrals
    double integral_surface_Ni_hex8(int i, int face_id, double x_hex8_arr[8], double y_hex8_arr[8], double z_hex8_arr[8]);

    // constructor
    IntegralHex8Class()
    {

    }
    IntegralHex8Class(GridHex8Struct &gh8s_in)
    {
        gh8s = gh8s_in;
    }

};

void IntegralHex8Class::evaluate_test_functions_derivatives()
{

    // integration points
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[8] = {-M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[8] = {-M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3};
    double c_arr[8] = {+M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};

    // iterate for each element
    for (int m = 0; m < gh8s.num_element; m++)
    {

        // initialize
        std::vector<double> jacobian_determinant_part_ml_vec;
        std::vector<std::vector<double>> N_part_ml_vec;
        std::vector<std::vector<double>> derivative_N_x_part_ml_vec;
        std::vector<std::vector<double>> derivative_N_y_part_ml_vec;
        std::vector<std::vector<double>> derivative_N_z_part_ml_vec;

        // get id of points around element
        int n00 = gh8s.element_p00_id_vec[m];
        int n01 = gh8s.element_p01_id_vec[m];
        int n02 = gh8s.element_p02_id_vec[m];
        int n03 = gh8s.element_p03_id_vec[m];
        int n04 = gh8s.element_p04_id_vec[m];
        int n05 = gh8s.element_p05_id_vec[m];
        int n06 = gh8s.element_p06_id_vec[m];
        int n07 = gh8s.element_p07_id_vec[m];

        // unpack x values
        double x00 = gh8s.point_pos_x_vec[n00];
        double x01 = gh8s.point_pos_x_vec[n01];
        double x02 = gh8s.point_pos_x_vec[n02];
        double x03 = gh8s.point_pos_x_vec[n03];
        double x04 = gh8s.point_pos_x_vec[n04];
        double x05 = gh8s.point_pos_x_vec[n05];
        double x06 = gh8s.point_pos_x_vec[n06];
        double x07 = gh8s.point_pos_x_vec[n07];

        // unpack y values
        double y00 = gh8s.point_pos_y_vec[n00];
        double y01 = gh8s.point_pos_y_vec[n01];
        double y02 = gh8s.point_pos_y_vec[n02];
        double y03 = gh8s.point_pos_y_vec[n03];
        double y04 = gh8s.point_pos_y_vec[n04];
        double y05 = gh8s.point_pos_y_vec[n05];
        double y06 = gh8s.point_pos_y_vec[n06];
        double y07 = gh8s.point_pos_y_vec[n07];

        // unpack z values
        double z00 = gh8s.point_pos_z_vec[n00];
        double z01 = gh8s.point_pos_z_vec[n01];
        double z02 = gh8s.point_pos_z_vec[n02];
        double z03 = gh8s.point_pos_z_vec[n03];
        double z04 = gh8s.point_pos_z_vec[n04];
        double z05 = gh8s.point_pos_z_vec[n05];
        double z06 = gh8s.point_pos_z_vec[n06];
        double z07 = gh8s.point_pos_z_vec[n07];

        // iterate for each integration point
        for (int l = 0; l < 8; l++)
        {
            // initialize
            std::vector<double> N_part_mli_vec;
            std::vector<double> derivative_N_x_part_mli_vec;
            std::vector<double> derivative_N_y_part_mli_vec;
            std::vector<double> derivative_N_z_part_mli_vec;

            // get a, b, and c values where function is evaluated
            double a = a_arr[l];
            double b = b_arr[l];
            double c = c_arr[l];

            // get derivatives of x and y with respect to a, b, and c
            double derivative_x_a = -0.125*(b*c*x00 - b*c*x01 + b*c*x02 - b*c*x03 - b*c*x04 + b*c*x05 - b*c*x06 + b*c*x07 - b*x00 - b*x01 + b*x02 + b*x03 + b*x04 + b*x05 - b*x06 - b*x07 - c*x00 + c*x01 - c*x02 + c*x03 - c*x04 + c*x05 - c*x06 + c*x07 + x00 + x01 - x02 - x03 + x04 + x05 - x06 - x07);
            double derivative_x_b = -0.125*(a*c*x00 - a*c*x01 + a*c*x02 - a*c*x03 - a*c*x04 + a*c*x05 - a*c*x06 + a*c*x07 - a*x00 - a*x01 + a*x02 + a*x03 + a*x04 + a*x05 - a*x06 - a*x07 - c*x00 + c*x01 + c*x02 - c*x03 + c*x04 - c*x05 - c*x06 + c*x07 + x00 + x01 + x02 + x03 - x04 - x05 - x06 - x07);
            double derivative_x_c = -0.125*(a*b*x00 - a*b*x01 + a*b*x02 - a*b*x03 - a*b*x04 + a*b*x05 - a*b*x06 + a*b*x07 - a*x00 + a*x01 - a*x02 + a*x03 - a*x04 + a*x05 - a*x06 + a*x07 - b*x00 + b*x01 + b*x02 - b*x03 + b*x04 - b*x05 - b*x06 + b*x07 + x00 - x01 - x02 + x03 + x04 - x05 - x06 + x07);
            double derivative_y_a = -0.125*(b*c*y00 - b*c*y01 + b*c*y02 - b*c*y03 - b*c*y04 + b*c*y05 - b*c*y06 + b*c*y07 - b*y00 - b*y01 + b*y02 + b*y03 + b*y04 + b*y05 - b*y06 - b*y07 - c*y00 + c*y01 - c*y02 + c*y03 - c*y04 + c*y05 - c*y06 + c*y07 + y00 + y01 - y02 - y03 + y04 + y05 - y06 - y07);
            double derivative_y_b = -0.125*(a*c*y00 - a*c*y01 + a*c*y02 - a*c*y03 - a*c*y04 + a*c*y05 - a*c*y06 + a*c*y07 - a*y00 - a*y01 + a*y02 + a*y03 + a*y04 + a*y05 - a*y06 - a*y07 - c*y00 + c*y01 + c*y02 - c*y03 + c*y04 - c*y05 - c*y06 + c*y07 + y00 + y01 + y02 + y03 - y04 - y05 - y06 - y07);
            double derivative_y_c = -0.125*(a*b*y00 - a*b*y01 + a*b*y02 - a*b*y03 - a*b*y04 + a*b*y05 - a*b*y06 + a*b*y07 - a*y00 + a*y01 - a*y02 + a*y03 - a*y04 + a*y05 - a*y06 + a*y07 - b*y00 + b*y01 + b*y02 - b*y03 + b*y04 - b*y05 - b*y06 + b*y07 + y00 - y01 - y02 + y03 + y04 - y05 - y06 + y07);
            double derivative_z_a = -0.125*(b*c*z00 - b*c*z01 + b*c*z02 - b*c*z03 - b*c*z04 + b*c*z05 - b*c*z06 + b*c*z07 - b*z00 - b*z01 + b*z02 + b*z03 + b*z04 + b*z05 - b*z06 - b*z07 - c*z00 + c*z01 - c*z02 + c*z03 - c*z04 + c*z05 - c*z06 + c*z07 + z00 + z01 - z02 - z03 + z04 + z05 - z06 - z07);
            double derivative_z_b = -0.125*(a*c*z00 - a*c*z01 + a*c*z02 - a*c*z03 - a*c*z04 + a*c*z05 - a*c*z06 + a*c*z07 - a*z00 - a*z01 + a*z02 + a*z03 + a*z04 + a*z05 - a*z06 - a*z07 - c*z00 + c*z01 + c*z02 - c*z03 + c*z04 - c*z05 - c*z06 + c*z07 + z00 + z01 + z02 + z03 - z04 - z05 - z06 - z07);
            double derivative_z_c = -0.125*(a*b*z00 - a*b*z01 + a*b*z02 - a*b*z03 - a*b*z04 + a*b*z05 - a*b*z06 + a*b*z07 - a*z00 + a*z01 - a*z02 + a*z03 - a*z04 + a*z05 - a*z06 + a*z07 - b*z00 + b*z01 + b*z02 - b*z03 + b*z04 - b*z05 - b*z06 + b*z07 + z00 - z01 - z02 + z03 + z04 - z05 - z06 + z07);

            // get jacobian and its inverse and determinant
            Eigen::Matrix3d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_x_c, derivative_y_a, derivative_y_b, derivative_y_c, derivative_z_a, derivative_z_b, derivative_z_c;
            Eigen::Matrix3d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int i = 0; i < 8; i++)
            {
        
                // get test function N
                double N = 0.;
                switch (i)
                {
                    case 0: N = 0.125*(1 - a)*(1 - b)*(1 - c); break;
                    case 1: N = 0.125*(1 - a)*(1 - b)*(1 + c); break;
                    case 2: N = 0.125*(1 + a)*(1 - b)*(1 + c); break;
                    case 3: N = 0.125*(1 + a)*(1 - b)*(1 - c); break;
                    case 4: N = 0.125*(1 - a)*(1 + b)*(1 - c); break;
                    case 5: N = 0.125*(1 - a)*(1 + b)*(1 + c); break;
                    case 6: N = 0.125*(1 + a)*(1 + b)*(1 + c); break;
                    case 7: N = 0.125*(1 + a)*(1 + b)*(1 - c); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                double derivative_N_b = 0.;
                double derivative_N_c = 0.;
                switch (i)
                {
                    case 0: derivative_N_a = -0.125*(b - 1)*(c - 1); derivative_N_b = -0.125*(a - 1)*(c - 1); derivative_N_c = -0.125*(a - 1)*(b - 1); break;
                    case 1: derivative_N_a = +0.125*(b - 1)*(c + 1); derivative_N_b = +0.125*(a - 1)*(c + 1); derivative_N_c = +0.125*(a - 1)*(b - 1); break;
                    case 2: derivative_N_a = -0.125*(b - 1)*(c + 1); derivative_N_b = -0.125*(a + 1)*(c + 1); derivative_N_c = -0.125*(a + 1)*(b - 1); break;
                    case 3: derivative_N_a = +0.125*(b - 1)*(c - 1); derivative_N_b = +0.125*(a + 1)*(c - 1); derivative_N_c = +0.125*(a + 1)*(b - 1); break;
                    case 4: derivative_N_a = +0.125*(b + 1)*(c - 1); derivative_N_b = +0.125*(a - 1)*(c - 1); derivative_N_c = +0.125*(a - 1)*(b + 1); break;
                    case 5: derivative_N_a = -0.125*(b + 1)*(c + 1); derivative_N_b = -0.125*(a - 1)*(c + 1); derivative_N_c = -0.125*(a - 1)*(b + 1); break;
                    case 6: derivative_N_a = +0.125*(b + 1)*(c + 1); derivative_N_b = +0.125*(a + 1)*(c + 1); derivative_N_c = +0.125*(a + 1)*(b + 1); break;
                    case 7: derivative_N_a = -0.125*(b + 1)*(c - 1); derivative_N_b = -0.125*(a + 1)*(c - 1); derivative_N_c = -0.125*(a + 1)*(b + 1); break;
                }

                // get vectors with derivatives of test functions wrt a, b, and c
                Eigen::RowVector3d derivative_N_abc_vec;
                derivative_N_abc_vec << derivative_N_a, derivative_N_b, derivative_N_c;

                // get vectors with derivatives of test functions wrt x, y, and z
                Eigen::RowVector3d derivative_N_xyz_vec;
                derivative_N_xyz_vec = derivative_N_abc_vec*jacobian_inverse_mat;
                double derivative_N_x = derivative_N_xyz_vec.coeffRef(0);
                double derivative_N_y = derivative_N_xyz_vec.coeffRef(1);
                double derivative_N_z = derivative_N_xyz_vec.coeffRef(2);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_N_x_part_mli_vec.push_back(derivative_N_x);
                derivative_N_y_part_mli_vec.push_back(derivative_N_y);
                derivative_N_z_part_mli_vec.push_back(derivative_N_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);
            derivative_N_y_part_ml_vec.push_back(derivative_N_y_part_mli_vec);
            derivative_N_z_part_ml_vec.push_back(derivative_N_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        N_vec.push_back(N_part_ml_vec);
        derivative_N_x_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_N_y_vec.push_back(derivative_N_y_part_ml_vec);
        derivative_N_z_vec.push_back(derivative_N_z_part_ml_vec);
        
    }

}

void IntegralHex8Class::evaluate_integral_div_Ni_hex8_dot_div_Nj_hex8()
{
    
    // iterate for each element
    for (int m = 0; m < gh8s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 8; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 8; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * (derivative_N_x_vec[m][l][i]*derivative_N_x_vec[m][l][j] + derivative_N_y_vec[m][l][i]*derivative_N_y_vec[m][l][j] + derivative_N_z_vec[m][l][i]*derivative_N_z_vec[m][l][j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_hex8_dot_div_Nj_hex8_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Class::evaluate_integral_Ni_hex8()
{
    
    // iterate for each element
    for (int m = 0; m < gh8s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<double> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 8; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_hex8_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Class::evaluate_integral_Ni_hex8_Nj_hex8()
{
    
    // iterate for each element
    for (int m = 0; m < gh8s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 8; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 8; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * N_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_hex8_Nj_hex8_vec.push_back(integral_part_i_vec);

    }

}

double IntegralHex8Class::integral_surface_Ni_hex8(int i, int face_id, double x_hex8_arr[8], double y_hex8_arr[8], double z_hex8_arr[8])
{

    // unpack x values
    double x00 = x_hex8_arr[0];
    double x01 = x_hex8_arr[1];
    double x02 = x_hex8_arr[2];
    double x03 = x_hex8_arr[3];
    double x04 = x_hex8_arr[4];
    double x05 = x_hex8_arr[5];
    double x06 = x_hex8_arr[6];
    double x07 = x_hex8_arr[7];

    // unpack y values
    double y00 = y_hex8_arr[0];
    double y01 = y_hex8_arr[1];
    double y02 = y_hex8_arr[2];
    double y03 = y_hex8_arr[3];
    double y04 = y_hex8_arr[4];
    double y05 = y_hex8_arr[5];
    double y06 = y_hex8_arr[6];
    double y07 = y_hex8_arr[7];

    // unpack z values
    double z00 = z_hex8_arr[0];
    double z01 = z_hex8_arr[1];
    double z02 = z_hex8_arr[2];
    double z03 = z_hex8_arr[3];
    double z04 = z_hex8_arr[4];
    double z05 = z_hex8_arr[5];
    double z06 = z_hex8_arr[6];
    double z07 = z_hex8_arr[7];

    // initialize integration points
    const double M_1_SQRT_3 = 1./sqrt(3);
    std::vector<double> sample_x_vec = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    std::vector<double> sample_y_vec = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};
    std::vector<double> sample_positive_vec = {+1, +1, +1, +1};
    std::vector<double> sample_negative_vec = {-1, -1, -1, -1};

    // get integration points depending on face
    std::vector<double> a_vec;
    std::vector<double> b_vec;
    std::vector<double> c_vec;
    switch (face_id)
    {
        case 0: a_vec = sample_negative_vec; b_vec = sample_x_vec; c_vec = sample_y_vec; break;  // -a
        case 1: a_vec = sample_positive_vec; b_vec = sample_x_vec; c_vec = sample_y_vec; break;  // +a
        case 2: a_vec = sample_x_vec; b_vec = sample_negative_vec; c_vec = sample_y_vec; break;  // -b
        case 3: a_vec = sample_x_vec; b_vec = sample_positive_vec; c_vec = sample_y_vec; break;  // +b
        case 4: a_vec = sample_x_vec; b_vec = sample_y_vec; c_vec = sample_negative_vec; break;  // -c
        case 5: a_vec = sample_x_vec; b_vec = sample_y_vec; c_vec = sample_positive_vec; break;  // +c
    }

    // initialize integral value
    double integral_value = 0.;

    // calculate for each point
    for (int m = 0; m < 4; m++)
    {

        // get a, b, and c values where function is evaluated
        double a = a_vec[m];
        double b = b_vec[m];
        double c = c_vec[m];

        // get test function Ni
        double Ni = 0.;
        switch (i)
        {
            case 0: Ni = 0.125*(1 - a)*(1 - b)*(1 - c); break;
            case 1: Ni = 0.125*(1 - a)*(1 - b)*(1 + c); break;
            case 2: Ni = 0.125*(1 + a)*(1 - b)*(1 + c); break;
            case 3: Ni = 0.125*(1 + a)*(1 - b)*(1 - c); break;
            case 4: Ni = 0.125*(1 - a)*(1 + b)*(1 - c); break;
            case 5: Ni = 0.125*(1 - a)*(1 + b)*(1 + c); break;
            case 6: Ni = 0.125*(1 + a)*(1 + b)*(1 + c); break;
            case 7: Ni = 0.125*(1 + a)*(1 + b)*(1 - c); break;
        }

        // get derivatives of x and y with respect to a, b, and c
        double derivative_x_a = -0.125*(b*c*x00 - b*c*x01 + b*c*x02 - b*c*x03 - b*c*x04 + b*c*x05 - b*c*x06 + b*c*x07 - b*x00 - b*x01 + b*x02 + b*x03 + b*x04 + b*x05 - b*x06 - b*x07 - c*x00 + c*x01 - c*x02 + c*x03 - c*x04 + c*x05 - c*x06 + c*x07 + x00 + x01 - x02 - x03 + x04 + x05 - x06 - x07);
        double derivative_x_b = -0.125*(a*c*x00 - a*c*x01 + a*c*x02 - a*c*x03 - a*c*x04 + a*c*x05 - a*c*x06 + a*c*x07 - a*x00 - a*x01 + a*x02 + a*x03 + a*x04 + a*x05 - a*x06 - a*x07 - c*x00 + c*x01 + c*x02 - c*x03 + c*x04 - c*x05 - c*x06 + c*x07 + x00 + x01 + x02 + x03 - x04 - x05 - x06 - x07);
        double derivative_x_c = -0.125*(a*b*x00 - a*b*x01 + a*b*x02 - a*b*x03 - a*b*x04 + a*b*x05 - a*b*x06 + a*b*x07 - a*x00 + a*x01 - a*x02 + a*x03 - a*x04 + a*x05 - a*x06 + a*x07 - b*x00 + b*x01 + b*x02 - b*x03 + b*x04 - b*x05 - b*x06 + b*x07 + x00 - x01 - x02 + x03 + x04 - x05 - x06 + x07);
        double derivative_y_a = -0.125*(b*c*y00 - b*c*y01 + b*c*y02 - b*c*y03 - b*c*y04 + b*c*y05 - b*c*y06 + b*c*y07 - b*y00 - b*y01 + b*y02 + b*y03 + b*y04 + b*y05 - b*y06 - b*y07 - c*y00 + c*y01 - c*y02 + c*y03 - c*y04 + c*y05 - c*y06 + c*y07 + y00 + y01 - y02 - y03 + y04 + y05 - y06 - y07);
        double derivative_y_b = -0.125*(a*c*y00 - a*c*y01 + a*c*y02 - a*c*y03 - a*c*y04 + a*c*y05 - a*c*y06 + a*c*y07 - a*y00 - a*y01 + a*y02 + a*y03 + a*y04 + a*y05 - a*y06 - a*y07 - c*y00 + c*y01 + c*y02 - c*y03 + c*y04 - c*y05 - c*y06 + c*y07 + y00 + y01 + y02 + y03 - y04 - y05 - y06 - y07);
        double derivative_y_c = -0.125*(a*b*y00 - a*b*y01 + a*b*y02 - a*b*y03 - a*b*y04 + a*b*y05 - a*b*y06 + a*b*y07 - a*y00 + a*y01 - a*y02 + a*y03 - a*y04 + a*y05 - a*y06 + a*y07 - b*y00 + b*y01 + b*y02 - b*y03 + b*y04 - b*y05 - b*y06 + b*y07 + y00 - y01 - y02 + y03 + y04 - y05 - y06 + y07);
        double derivative_z_a = -0.125*(b*c*z00 - b*c*z01 + b*c*z02 - b*c*z03 - b*c*z04 + b*c*z05 - b*c*z06 + b*c*z07 - b*z00 - b*z01 + b*z02 + b*z03 + b*z04 + b*z05 - b*z06 - b*z07 - c*z00 + c*z01 - c*z02 + c*z03 - c*z04 + c*z05 - c*z06 + c*z07 + z00 + z01 - z02 - z03 + z04 + z05 - z06 - z07);
        double derivative_z_b = -0.125*(a*c*z00 - a*c*z01 + a*c*z02 - a*c*z03 - a*c*z04 + a*c*z05 - a*c*z06 + a*c*z07 - a*z00 - a*z01 + a*z02 + a*z03 + a*z04 + a*z05 - a*z06 - a*z07 - c*z00 + c*z01 + c*z02 - c*z03 + c*z04 - c*z05 - c*z06 + c*z07 + z00 + z01 + z02 + z03 - z04 - z05 - z06 - z07);
        double derivative_z_c = -0.125*(a*b*z00 - a*b*z01 + a*b*z02 - a*b*z03 - a*b*z04 + a*b*z05 - a*b*z06 + a*b*z07 - a*z00 + a*z01 - a*z02 + a*z03 - a*z04 + a*z05 - a*z06 + a*z07 - b*z00 + b*z01 + b*z02 - b*z03 + b*z04 - b*z05 - b*z06 + b*z07 + z00 - z01 - z02 + z03 + z04 - z05 - z06 + z07);

        // get jacobian and its inverse and determinant
        // jacobian_mat << derivative_x_a, derivative_x_b, derivative_x_c,
        Eigen::Matrix2d jacobian_mat;
        switch (face_id)
        {
            case 0: jacobian_mat << derivative_y_b, derivative_y_c, derivative_z_b, derivative_z_c; break;  // -a
            case 1: jacobian_mat << derivative_y_b, derivative_y_c, derivative_z_b, derivative_z_c; break;  // +a
            case 2: jacobian_mat << derivative_x_a, derivative_x_c, derivative_z_a, derivative_z_c; break;  // -b
            case 3: jacobian_mat << derivative_x_a, derivative_x_c, derivative_z_a, derivative_z_c; break;  // +b
            case 4: jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b; break;  // -c
            case 5: jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b; break;  // +c
        }
        double jacobian_determinant = jacobian_mat.determinant();

        // evaluate part of integral value
        integral_value += jacobian_determinant * Ni;

    }

    return integral_value;

}

#endif
