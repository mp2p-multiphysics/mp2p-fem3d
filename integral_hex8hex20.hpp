#ifndef INTEGRAL_HEX8HEX20
#define INTEGRAL_HEX8HEX20
#include <vector>
#include "Eigen/Eigen"

class IntegralHex8Hex20Class
{

    public:

    // variables
    GridHex20Struct gh20s;
    GridHex8Struct gh8s;

    // vectors with test functions and derivatives
    std::vector<std::vector<double>> jacobian_determinant_vec;
    std::vector<std::vector<std::vector<double>>> M_vec;
    std::vector<std::vector<std::vector<double>>> derivative_M_x_vec;
    std::vector<std::vector<std::vector<double>>> derivative_M_y_vec;
    std::vector<std::vector<std::vector<double>>> derivative_M_z_vec;
    std::vector<std::vector<std::vector<double>>> N_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_x_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_y_vec;
    std::vector<std::vector<std::vector<double>>> derivative_N_z_vec;

    // vectors with integrals
    std::vector<std::vector<double>> integral_Mi_hex20_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_hex8_derivative_Mj_hex20_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_hex8_derivative_Mj_hex20_y_vec;
    std::vector<std::vector<std::vector<double>>> integral_Ni_hex8_derivative_Mj_hex20_z_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_hex20_derivative_Nj_hex8_x_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_hex20_derivative_Nj_hex8_y_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_hex20_derivative_Nj_hex8_z_vec;
    std::vector<std::vector<std::vector<double>>> integral_div_Mi_hex20_dot_div_Mj_hex20_vec;
    std::vector<std::vector<std::vector<double>>> integral_Mi_hex20_Mj_hex20_vec;
    std::vector<std::vector<std::vector<std::vector<double>>>> integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x_vec;
    std::vector<std::vector<std::vector<std::vector<double>>>> integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y_vec;
    std::vector<std::vector<std::vector<std::vector<double>>>> integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z_vec;

    // functions for computing integrals
    void evaluate_test_functions_derivatives();
    void evaluate_integral_Mi_hex20();
    void evaluate_integral_Ni_hex8_derivative_Mj_hex20_x();
    void evaluate_integral_Ni_hex8_derivative_Mj_hex20_y();
    void evaluate_integral_Ni_hex8_derivative_Mj_hex20_z();
    void evaluate_integral_Mi_hex20_derivative_Nj_hex8_x();
    void evaluate_integral_Mi_hex20_derivative_Nj_hex8_y();
    void evaluate_integral_Mi_hex20_derivative_Nj_hex8_z();
    void evaluate_integral_div_Mi_hex20_dot_div_Mj_hex20();
    void evaluate_integral_Mi_hex20_Mj_hex20();
    void evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x();
    void evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y();
    void evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z();

    // constructor
    IntegralHex8Hex20Class()
    {

    }
    IntegralHex8Hex20Class(GridHex20Struct &gh20s_in, GridHex8Struct &gh8s_in)
    {
        gh20s = gh20s_in;
        gh8s = gh8s_in;
    }

};

void IntegralHex8Hex20Class::evaluate_test_functions_derivatives()
{

    // integration points
    const double M_SQRT_3_5 = sqrt(0.6);
    const double M_125_729 = 125./729.;
    const double M_200_729 = 200./729.;
    const double M_320_729 = 320./729.;
    const double M_512_729 = 512./729.;
    double a_arr[27] = {
        -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5,           0,
        -M_SQRT_3_5, +M_SQRT_3_5,           0, -M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5,           0, -M_SQRT_3_5,
        +M_SQRT_3_5,           0,           0,           0, -M_SQRT_3_5, +M_SQRT_3_5,           0,           0,           0
    };
    double b_arr[27] = {
        -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5,
                  0,           0, +M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5,           0,
                  0, +M_SQRT_3_5,           0, -M_SQRT_3_5,           0,           0, +M_SQRT_3_5,           0,           0
    };
    double c_arr[27] = {
        -M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5,
        -M_SQRT_3_5, -M_SQRT_3_5, -M_SQRT_3_5,           0,           0,           0,           0, +M_SQRT_3_5, +M_SQRT_3_5,
        +M_SQRT_3_5, +M_SQRT_3_5, -M_SQRT_3_5,           0,           0,           0,           0, +M_SQRT_3_5,           0
    };
    double w_arr[27] = {
        M_125_729, M_125_729, M_125_729, M_125_729, M_125_729, M_125_729, M_125_729, M_125_729, M_200_729,
        M_200_729, M_200_729, M_200_729, M_200_729, M_200_729, M_200_729, M_200_729, M_200_729, M_200_729,
        M_200_729, M_200_729, M_320_729, M_320_729, M_320_729, M_320_729, M_320_729, M_320_729, M_512_729
    };

    // iterate for each hex20 element
    for (int m = 0; m < gh20s.num_element; m++)
    {

        // initialize
        std::vector<double> jacobian_determinant_part_ml_vec;
        std::vector<std::vector<double>> M_part_ml_vec;
        std::vector<std::vector<double>> derivative_M_x_part_ml_vec;
        std::vector<std::vector<double>> derivative_M_y_part_ml_vec;
        std::vector<std::vector<double>> derivative_M_z_part_ml_vec;

        // get id of points around hex20 element
        double n00 = gh20s.element_p00_id_vec[m]; double n10 = gh20s.element_p10_id_vec[m];
        double n01 = gh20s.element_p01_id_vec[m]; double n11 = gh20s.element_p11_id_vec[m];
        double n02 = gh20s.element_p02_id_vec[m]; double n12 = gh20s.element_p12_id_vec[m];
        double n03 = gh20s.element_p03_id_vec[m]; double n13 = gh20s.element_p13_id_vec[m];
        double n04 = gh20s.element_p04_id_vec[m]; double n14 = gh20s.element_p14_id_vec[m];
        double n05 = gh20s.element_p05_id_vec[m]; double n15 = gh20s.element_p15_id_vec[m];
        double n06 = gh20s.element_p06_id_vec[m]; double n16 = gh20s.element_p16_id_vec[m];
        double n07 = gh20s.element_p07_id_vec[m]; double n17 = gh20s.element_p17_id_vec[m];
        double n08 = gh20s.element_p08_id_vec[m]; double n18 = gh20s.element_p18_id_vec[m];
        double n09 = gh20s.element_p09_id_vec[m]; double n19 = gh20s.element_p19_id_vec[m];

        // unpack hex20 x values
        double x00 = gh20s.point_pos_x_vec[n00]; double x10 = gh20s.point_pos_x_vec[n10];
        double x01 = gh20s.point_pos_x_vec[n01]; double x11 = gh20s.point_pos_x_vec[n11];
        double x02 = gh20s.point_pos_x_vec[n02]; double x12 = gh20s.point_pos_x_vec[n12];
        double x03 = gh20s.point_pos_x_vec[n03]; double x13 = gh20s.point_pos_x_vec[n13];
        double x04 = gh20s.point_pos_x_vec[n04]; double x14 = gh20s.point_pos_x_vec[n14];
        double x05 = gh20s.point_pos_x_vec[n05]; double x15 = gh20s.point_pos_x_vec[n15];
        double x06 = gh20s.point_pos_x_vec[n06]; double x16 = gh20s.point_pos_x_vec[n16];
        double x07 = gh20s.point_pos_x_vec[n07]; double x17 = gh20s.point_pos_x_vec[n17];
        double x08 = gh20s.point_pos_x_vec[n08]; double x18 = gh20s.point_pos_x_vec[n18];
        double x09 = gh20s.point_pos_x_vec[n09]; double x19 = gh20s.point_pos_x_vec[n19];

        // unpack hex20 y values
        double y00 = gh20s.point_pos_y_vec[n00]; double y10 = gh20s.point_pos_y_vec[n10];
        double y01 = gh20s.point_pos_y_vec[n01]; double y11 = gh20s.point_pos_y_vec[n11];
        double y02 = gh20s.point_pos_y_vec[n02]; double y12 = gh20s.point_pos_y_vec[n12];
        double y03 = gh20s.point_pos_y_vec[n03]; double y13 = gh20s.point_pos_y_vec[n13];
        double y04 = gh20s.point_pos_y_vec[n04]; double y14 = gh20s.point_pos_y_vec[n14];
        double y05 = gh20s.point_pos_y_vec[n05]; double y15 = gh20s.point_pos_y_vec[n15];
        double y06 = gh20s.point_pos_y_vec[n06]; double y16 = gh20s.point_pos_y_vec[n16];
        double y07 = gh20s.point_pos_y_vec[n07]; double y17 = gh20s.point_pos_y_vec[n17];
        double y08 = gh20s.point_pos_y_vec[n08]; double y18 = gh20s.point_pos_y_vec[n18];
        double y09 = gh20s.point_pos_y_vec[n09]; double y19 = gh20s.point_pos_y_vec[n19];

        // unpack hex20 z values
        double z00 = gh20s.point_pos_z_vec[n00]; double z10 = gh20s.point_pos_z_vec[n10];
        double z01 = gh20s.point_pos_z_vec[n01]; double z11 = gh20s.point_pos_z_vec[n11];
        double z02 = gh20s.point_pos_z_vec[n02]; double z12 = gh20s.point_pos_z_vec[n12];
        double z03 = gh20s.point_pos_z_vec[n03]; double z13 = gh20s.point_pos_z_vec[n13];
        double z04 = gh20s.point_pos_z_vec[n04]; double z14 = gh20s.point_pos_z_vec[n14];
        double z05 = gh20s.point_pos_z_vec[n05]; double z15 = gh20s.point_pos_z_vec[n15];
        double z06 = gh20s.point_pos_z_vec[n06]; double z16 = gh20s.point_pos_z_vec[n16];
        double z07 = gh20s.point_pos_z_vec[n07]; double z17 = gh20s.point_pos_z_vec[n17];
        double z08 = gh20s.point_pos_z_vec[n08]; double z18 = gh20s.point_pos_z_vec[n18];
        double z09 = gh20s.point_pos_z_vec[n09]; double z19 = gh20s.point_pos_z_vec[n19];

        // iterate for each integration point
        for (int l = 0; l < 27; l++)
        {

            // initialize
            std::vector<double> M_part_mli_vec;
            std::vector<double> derivative_M_x_part_mli_vec;
            std::vector<double> derivative_M_y_part_mli_vec;
            std::vector<double> derivative_M_z_part_mli_vec;

            // get a, b, and c values where function is evaluated
            double a = a_arr[l];
            double b = b_arr[l];
            double c = c_arr[l];

            // get derivatives of x and y with respect to a, b, and c
            double derivative_x_a = 0.125*(2*a*b*c*x00 - 2*a*b*c*x02 + 4*a*b*c*x03 - 2*a*b*c*x04 + 2*a*b*c*x06 - 4*a*b*c*x07 - 2*a*b*c*x12 + 2*a*b*c*x14 - 4*a*b*c*x15 + 2*a*b*c*x16 - 2*a*b*c*x18 + 4*a*b*c*x19 - 2*a*b*x00 - 2*a*b*x02 + 4*a*b*x03 - 2*a*b*x04 - 2*a*b*x06 + 4*a*b*x07 + 2*a*b*x12 + 2*a*b*x14 - 4*a*b*x15 + 2*a*b*x16 + 2*a*b*x18 - 4*a*b*x19 - 2*a*c*x00 + 2*a*c*x02 - 4*a*c*x03 + 2*a*c*x04 - 2*a*c*x06 + 4*a*c*x07 - 2*a*c*x12 + 2*a*c*x14 - 4*a*c*x15 + 2*a*c*x16 - 2*a*c*x18 + 4*a*c*x19 + 2*a*x00 + 2*a*x02 - 4*a*x03 + 2*a*x04 + 2*a*x06 - 4*a*x07 + 2*a*x12 + 2*a*x14 - 4*a*x15 + 2*a*x16 + 2*a*x18 - 4*a*x19 + b*b*c*x00 - b*b*c*x02 + b*b*c*x04 - b*b*c*x06 - 2*b*b*c*x08 + 2*b*b*c*x09 - 2*b*b*c*x10 + 2*b*b*c*x11 + b*b*c*x12 - b*b*c*x14 + b*b*c*x16 - b*b*c*x18 - b*b*x00 - b*b*x02 + b*b*x04 + b*b*x06 + 2*b*b*x08 + 2*b*b*x09 - 2*b*b*x10 - 2*b*b*x11 - b*b*x12 - b*b*x14 + b*b*x16 + b*b*x18 + b*c*c*x00 - 2*b*c*c*x01 + b*c*c*x02 - b*c*c*x04 + 2*b*c*c*x05 - b*c*c*x06 - b*c*c*x12 + 2*b*c*c*x13 - b*c*c*x14 + b*c*c*x16 - 2*b*c*c*x17 + b*c*c*x18 - b*c*x00 + b*c*x02 - b*c*x04 + b*c*x06 + b*c*x12 - b*c*x14 + b*c*x16 - b*c*x18 + 2*b*x01 - 2*b*x05 - 2*b*x13 + 2*b*x17 - c*c*x00 + 2*c*c*x01 - c*c*x02 + c*c*x04 - 2*c*c*x05 + c*c*x06 - c*c*x12 + 2*c*c*x13 - c*c*x14 + c*c*x16 - 2*c*c*x17 + c*c*x18 + 2*c*x08 - 2*c*x09 + 2*c*x10 - 2*c*x11 + x00 - 2*x01 + x02 - x04 + 2*x05 - x06 - 2*x08 - 2*x09 + 2*x10 + 2*x11 + x12 - 2*x13 + x14 - x16 + 2*x17 - x18);
            double derivative_x_b = 0.125*(a*a*c*x00 - a*a*c*x02 + 2*a*a*c*x03 - a*a*c*x04 + a*a*c*x06 - 2*a*a*c*x07 - a*a*c*x12 + a*a*c*x14 - 2*a*a*c*x15 + a*a*c*x16 - a*a*c*x18 + 2*a*a*c*x19 - a*a*x00 - a*a*x02 + 2*a*a*x03 - a*a*x04 - a*a*x06 + 2*a*a*x07 + a*a*x12 + a*a*x14 - 2*a*a*x15 + a*a*x16 + a*a*x18 - 2*a*a*x19 + 2*a*b*c*x00 - 2*a*b*c*x02 + 2*a*b*c*x04 - 2*a*b*c*x06 - 4*a*b*c*x08 + 4*a*b*c*x09 - 4*a*b*c*x10 + 4*a*b*c*x11 + 2*a*b*c*x12 - 2*a*b*c*x14 + 2*a*b*c*x16 - 2*a*b*c*x18 - 2*a*b*x00 - 2*a*b*x02 + 2*a*b*x04 + 2*a*b*x06 + 4*a*b*x08 + 4*a*b*x09 - 4*a*b*x10 - 4*a*b*x11 - 2*a*b*x12 - 2*a*b*x14 + 2*a*b*x16 + 2*a*b*x18 + a*c*c*x00 - 2*a*c*c*x01 + a*c*c*x02 - a*c*c*x04 + 2*a*c*c*x05 - a*c*c*x06 - a*c*c*x12 + 2*a*c*c*x13 - a*c*c*x14 + a*c*c*x16 - 2*a*c*c*x17 + a*c*c*x18 - a*c*x00 + a*c*x02 - a*c*x04 + a*c*x06 + a*c*x12 - a*c*x14 + a*c*x16 - a*c*x18 + 2*a*x01 - 2*a*x05 - 2*a*x13 + 2*a*x17 - 2*b*c*x00 + 2*b*c*x02 + 2*b*c*x04 - 2*b*c*x06 + 4*b*c*x08 - 4*b*c*x09 - 4*b*c*x10 + 4*b*c*x11 - 2*b*c*x12 + 2*b*c*x14 + 2*b*c*x16 - 2*b*c*x18 + 2*b*x00 + 2*b*x02 + 2*b*x04 + 2*b*x06 - 4*b*x08 - 4*b*x09 - 4*b*x10 - 4*b*x11 + 2*b*x12 + 2*b*x14 + 2*b*x16 + 2*b*x18 - c*c*x00 + 2*c*c*x01 - c*c*x02 - c*c*x04 + 2*c*c*x05 - c*c*x06 + c*c*x12 - 2*c*c*x13 + c*c*x14 + c*c*x16 - 2*c*c*x17 + c*c*x18 - 2*c*x03 + 2*c*x07 + 2*c*x15 - 2*c*x19 + x00 - 2*x01 + x02 - 2*x03 + x04 - 2*x05 + x06 - 2*x07 - x12 + 2*x13 - x14 + 2*x15 - x16 + 2*x17 - x18 + 2*x19);
            double derivative_x_c = 0.125*(a*a*b*x00 - a*a*b*x02 + 2*a*a*b*x03 - a*a*b*x04 + a*a*b*x06 - 2*a*a*b*x07 - a*a*b*x12 + a*a*b*x14 - 2*a*a*b*x15 + a*a*b*x16 - a*a*b*x18 + 2*a*a*b*x19 - a*a*x00 + a*a*x02 - 2*a*a*x03 + a*a*x04 - a*a*x06 + 2*a*a*x07 - a*a*x12 + a*a*x14 - 2*a*a*x15 + a*a*x16 - a*a*x18 + 2*a*a*x19 + a*b*b*x00 - a*b*b*x02 + a*b*b*x04 - a*b*b*x06 - 2*a*b*b*x08 + 2*a*b*b*x09 - 2*a*b*b*x10 + 2*a*b*b*x11 + a*b*b*x12 - a*b*b*x14 + a*b*b*x16 - a*b*b*x18 + 2*a*b*c*x00 - 4*a*b*c*x01 + 2*a*b*c*x02 - 2*a*b*c*x04 + 4*a*b*c*x05 - 2*a*b*c*x06 - 2*a*b*c*x12 + 4*a*b*c*x13 - 2*a*b*c*x14 + 2*a*b*c*x16 - 4*a*b*c*x17 + 2*a*b*c*x18 - a*b*x00 + a*b*x02 - a*b*x04 + a*b*x06 + a*b*x12 - a*b*x14 + a*b*x16 - a*b*x18 - 2*a*c*x00 + 4*a*c*x01 - 2*a*c*x02 + 2*a*c*x04 - 4*a*c*x05 + 2*a*c*x06 - 2*a*c*x12 + 4*a*c*x13 - 2*a*c*x14 + 2*a*c*x16 - 4*a*c*x17 + 2*a*c*x18 + 2*a*x08 - 2*a*x09 + 2*a*x10 - 2*a*x11 - b*b*x00 + b*b*x02 + b*b*x04 - b*b*x06 + 2*b*b*x08 - 2*b*b*x09 - 2*b*b*x10 + 2*b*b*x11 - b*b*x12 + b*b*x14 + b*b*x16 - b*b*x18 - 2*b*c*x00 + 4*b*c*x01 - 2*b*c*x02 - 2*b*c*x04 + 4*b*c*x05 - 2*b*c*x06 + 2*b*c*x12 - 4*b*c*x13 + 2*b*c*x14 + 2*b*c*x16 - 4*b*c*x17 + 2*b*c*x18 - 2*b*x03 + 2*b*x07 + 2*b*x15 - 2*b*x19 + 2*c*x00 - 4*c*x01 + 2*c*x02 + 2*c*x04 - 4*c*x05 + 2*c*x06 + 2*c*x12 - 4*c*x13 + 2*c*x14 + 2*c*x16 - 4*c*x17 + 2*c*x18 + x00 - x02 + 2*x03 - x04 + x06 - 2*x07 - 2*x08 + 2*x09 + 2*x10 - 2*x11 + x12 - x14 + 2*x15 - x16 + x18 - 2*x19);
            double derivative_y_a = 0.125*(2*a*b*c*y00 - 2*a*b*c*y02 + 4*a*b*c*y03 - 2*a*b*c*y04 + 2*a*b*c*y06 - 4*a*b*c*y07 - 2*a*b*c*y12 + 2*a*b*c*y14 - 4*a*b*c*y15 + 2*a*b*c*y16 - 2*a*b*c*y18 + 4*a*b*c*y19 - 2*a*b*y00 - 2*a*b*y02 + 4*a*b*y03 - 2*a*b*y04 - 2*a*b*y06 + 4*a*b*y07 + 2*a*b*y12 + 2*a*b*y14 - 4*a*b*y15 + 2*a*b*y16 + 2*a*b*y18 - 4*a*b*y19 - 2*a*c*y00 + 2*a*c*y02 - 4*a*c*y03 + 2*a*c*y04 - 2*a*c*y06 + 4*a*c*y07 - 2*a*c*y12 + 2*a*c*y14 - 4*a*c*y15 + 2*a*c*y16 - 2*a*c*y18 + 4*a*c*y19 + 2*a*y00 + 2*a*y02 - 4*a*y03 + 2*a*y04 + 2*a*y06 - 4*a*y07 + 2*a*y12 + 2*a*y14 - 4*a*y15 + 2*a*y16 + 2*a*y18 - 4*a*y19 + b*b*c*y00 - b*b*c*y02 + b*b*c*y04 - b*b*c*y06 - 2*b*b*c*y08 + 2*b*b*c*y09 - 2*b*b*c*y10 + 2*b*b*c*y11 + b*b*c*y12 - b*b*c*y14 + b*b*c*y16 - b*b*c*y18 - b*b*y00 - b*b*y02 + b*b*y04 + b*b*y06 + 2*b*b*y08 + 2*b*b*y09 - 2*b*b*y10 - 2*b*b*y11 - b*b*y12 - b*b*y14 + b*b*y16 + b*b*y18 + b*c*c*y00 - 2*b*c*c*y01 + b*c*c*y02 - b*c*c*y04 + 2*b*c*c*y05 - b*c*c*y06 - b*c*c*y12 + 2*b*c*c*y13 - b*c*c*y14 + b*c*c*y16 - 2*b*c*c*y17 + b*c*c*y18 - b*c*y00 + b*c*y02 - b*c*y04 + b*c*y06 + b*c*y12 - b*c*y14 + b*c*y16 - b*c*y18 + 2*b*y01 - 2*b*y05 - 2*b*y13 + 2*b*y17 - c*c*y00 + 2*c*c*y01 - c*c*y02 + c*c*y04 - 2*c*c*y05 + c*c*y06 - c*c*y12 + 2*c*c*y13 - c*c*y14 + c*c*y16 - 2*c*c*y17 + c*c*y18 + 2*c*y08 - 2*c*y09 + 2*c*y10 - 2*c*y11 + y00 - 2*y01 + y02 - y04 + 2*y05 - y06 - 2*y08 - 2*y09 + 2*y10 + 2*y11 + y12 - 2*y13 + y14 - y16 + 2*y17 - y18);
            double derivative_y_b = 0.125*(a*a*c*y00 - a*a*c*y02 + 2*a*a*c*y03 - a*a*c*y04 + a*a*c*y06 - 2*a*a*c*y07 - a*a*c*y12 + a*a*c*y14 - 2*a*a*c*y15 + a*a*c*y16 - a*a*c*y18 + 2*a*a*c*y19 - a*a*y00 - a*a*y02 + 2*a*a*y03 - a*a*y04 - a*a*y06 + 2*a*a*y07 + a*a*y12 + a*a*y14 - 2*a*a*y15 + a*a*y16 + a*a*y18 - 2*a*a*y19 + 2*a*b*c*y00 - 2*a*b*c*y02 + 2*a*b*c*y04 - 2*a*b*c*y06 - 4*a*b*c*y08 + 4*a*b*c*y09 - 4*a*b*c*y10 + 4*a*b*c*y11 + 2*a*b*c*y12 - 2*a*b*c*y14 + 2*a*b*c*y16 - 2*a*b*c*y18 - 2*a*b*y00 - 2*a*b*y02 + 2*a*b*y04 + 2*a*b*y06 + 4*a*b*y08 + 4*a*b*y09 - 4*a*b*y10 - 4*a*b*y11 - 2*a*b*y12 - 2*a*b*y14 + 2*a*b*y16 + 2*a*b*y18 + a*c*c*y00 - 2*a*c*c*y01 + a*c*c*y02 - a*c*c*y04 + 2*a*c*c*y05 - a*c*c*y06 - a*c*c*y12 + 2*a*c*c*y13 - a*c*c*y14 + a*c*c*y16 - 2*a*c*c*y17 + a*c*c*y18 - a*c*y00 + a*c*y02 - a*c*y04 + a*c*y06 + a*c*y12 - a*c*y14 + a*c*y16 - a*c*y18 + 2*a*y01 - 2*a*y05 - 2*a*y13 + 2*a*y17 - 2*b*c*y00 + 2*b*c*y02 + 2*b*c*y04 - 2*b*c*y06 + 4*b*c*y08 - 4*b*c*y09 - 4*b*c*y10 + 4*b*c*y11 - 2*b*c*y12 + 2*b*c*y14 + 2*b*c*y16 - 2*b*c*y18 + 2*b*y00 + 2*b*y02 + 2*b*y04 + 2*b*y06 - 4*b*y08 - 4*b*y09 - 4*b*y10 - 4*b*y11 + 2*b*y12 + 2*b*y14 + 2*b*y16 + 2*b*y18 - c*c*y00 + 2*c*c*y01 - c*c*y02 - c*c*y04 + 2*c*c*y05 - c*c*y06 + c*c*y12 - 2*c*c*y13 + c*c*y14 + c*c*y16 - 2*c*c*y17 + c*c*y18 - 2*c*y03 + 2*c*y07 + 2*c*y15 - 2*c*y19 + y00 - 2*y01 + y02 - 2*y03 + y04 - 2*y05 + y06 - 2*y07 - y12 + 2*y13 - y14 + 2*y15 - y16 + 2*y17 - y18 + 2*y19);
            double derivative_y_c = 0.125*(a*a*b*y00 - a*a*b*y02 + 2*a*a*b*y03 - a*a*b*y04 + a*a*b*y06 - 2*a*a*b*y07 - a*a*b*y12 + a*a*b*y14 - 2*a*a*b*y15 + a*a*b*y16 - a*a*b*y18 + 2*a*a*b*y19 - a*a*y00 + a*a*y02 - 2*a*a*y03 + a*a*y04 - a*a*y06 + 2*a*a*y07 - a*a*y12 + a*a*y14 - 2*a*a*y15 + a*a*y16 - a*a*y18 + 2*a*a*y19 + a*b*b*y00 - a*b*b*y02 + a*b*b*y04 - a*b*b*y06 - 2*a*b*b*y08 + 2*a*b*b*y09 - 2*a*b*b*y10 + 2*a*b*b*y11 + a*b*b*y12 - a*b*b*y14 + a*b*b*y16 - a*b*b*y18 + 2*a*b*c*y00 - 4*a*b*c*y01 + 2*a*b*c*y02 - 2*a*b*c*y04 + 4*a*b*c*y05 - 2*a*b*c*y06 - 2*a*b*c*y12 + 4*a*b*c*y13 - 2*a*b*c*y14 + 2*a*b*c*y16 - 4*a*b*c*y17 + 2*a*b*c*y18 - a*b*y00 + a*b*y02 - a*b*y04 + a*b*y06 + a*b*y12 - a*b*y14 + a*b*y16 - a*b*y18 - 2*a*c*y00 + 4*a*c*y01 - 2*a*c*y02 + 2*a*c*y04 - 4*a*c*y05 + 2*a*c*y06 - 2*a*c*y12 + 4*a*c*y13 - 2*a*c*y14 + 2*a*c*y16 - 4*a*c*y17 + 2*a*c*y18 + 2*a*y08 - 2*a*y09 + 2*a*y10 - 2*a*y11 - b*b*y00 + b*b*y02 + b*b*y04 - b*b*y06 + 2*b*b*y08 - 2*b*b*y09 - 2*b*b*y10 + 2*b*b*y11 - b*b*y12 + b*b*y14 + b*b*y16 - b*b*y18 - 2*b*c*y00 + 4*b*c*y01 - 2*b*c*y02 - 2*b*c*y04 + 4*b*c*y05 - 2*b*c*y06 + 2*b*c*y12 - 4*b*c*y13 + 2*b*c*y14 + 2*b*c*y16 - 4*b*c*y17 + 2*b*c*y18 - 2*b*y03 + 2*b*y07 + 2*b*y15 - 2*b*y19 + 2*c*y00 - 4*c*y01 + 2*c*y02 + 2*c*y04 - 4*c*y05 + 2*c*y06 + 2*c*y12 - 4*c*y13 + 2*c*y14 + 2*c*y16 - 4*c*y17 + 2*c*y18 + y00 - y02 + 2*y03 - y04 + y06 - 2*y07 - 2*y08 + 2*y09 + 2*y10 - 2*y11 + y12 - y14 + 2*y15 - y16 + y18 - 2*y19);
            double derivative_z_a = 0.125*(2*a*b*c*z00 - 2*a*b*c*z02 + 4*a*b*c*z03 - 2*a*b*c*z04 + 2*a*b*c*z06 - 4*a*b*c*z07 - 2*a*b*c*z12 + 2*a*b*c*z14 - 4*a*b*c*z15 + 2*a*b*c*z16 - 2*a*b*c*z18 + 4*a*b*c*z19 - 2*a*b*z00 - 2*a*b*z02 + 4*a*b*z03 - 2*a*b*z04 - 2*a*b*z06 + 4*a*b*z07 + 2*a*b*z12 + 2*a*b*z14 - 4*a*b*z15 + 2*a*b*z16 + 2*a*b*z18 - 4*a*b*z19 - 2*a*c*z00 + 2*a*c*z02 - 4*a*c*z03 + 2*a*c*z04 - 2*a*c*z06 + 4*a*c*z07 - 2*a*c*z12 + 2*a*c*z14 - 4*a*c*z15 + 2*a*c*z16 - 2*a*c*z18 + 4*a*c*z19 + 2*a*z00 + 2*a*z02 - 4*a*z03 + 2*a*z04 + 2*a*z06 - 4*a*z07 + 2*a*z12 + 2*a*z14 - 4*a*z15 + 2*a*z16 + 2*a*z18 - 4*a*z19 + b*b*c*z00 - b*b*c*z02 + b*b*c*z04 - b*b*c*z06 - 2*b*b*c*z08 + 2*b*b*c*z09 - 2*b*b*c*z10 + 2*b*b*c*z11 + b*b*c*z12 - b*b*c*z14 + b*b*c*z16 - b*b*c*z18 - b*b*z00 - b*b*z02 + b*b*z04 + b*b*z06 + 2*b*b*z08 + 2*b*b*z09 - 2*b*b*z10 - 2*b*b*z11 - b*b*z12 - b*b*z14 + b*b*z16 + b*b*z18 + b*c*c*z00 - 2*b*c*c*z01 + b*c*c*z02 - b*c*c*z04 + 2*b*c*c*z05 - b*c*c*z06 - b*c*c*z12 + 2*b*c*c*z13 - b*c*c*z14 + b*c*c*z16 - 2*b*c*c*z17 + b*c*c*z18 - b*c*z00 + b*c*z02 - b*c*z04 + b*c*z06 + b*c*z12 - b*c*z14 + b*c*z16 - b*c*z18 + 2*b*z01 - 2*b*z05 - 2*b*z13 + 2*b*z17 - c*c*z00 + 2*c*c*z01 - c*c*z02 + c*c*z04 - 2*c*c*z05 + c*c*z06 - c*c*z12 + 2*c*c*z13 - c*c*z14 + c*c*z16 - 2*c*c*z17 + c*c*z18 + 2*c*z08 - 2*c*z09 + 2*c*z10 - 2*c*z11 + z00 - 2*z01 + z02 - z04 + 2*z05 - z06 - 2*z08 - 2*z09 + 2*z10 + 2*z11 + z12 - 2*z13 + z14 - z16 + 2*z17 - z18);
            double derivative_z_b = 0.125*(a*a*c*z00 - a*a*c*z02 + 2*a*a*c*z03 - a*a*c*z04 + a*a*c*z06 - 2*a*a*c*z07 - a*a*c*z12 + a*a*c*z14 - 2*a*a*c*z15 + a*a*c*z16 - a*a*c*z18 + 2*a*a*c*z19 - a*a*z00 - a*a*z02 + 2*a*a*z03 - a*a*z04 - a*a*z06 + 2*a*a*z07 + a*a*z12 + a*a*z14 - 2*a*a*z15 + a*a*z16 + a*a*z18 - 2*a*a*z19 + 2*a*b*c*z00 - 2*a*b*c*z02 + 2*a*b*c*z04 - 2*a*b*c*z06 - 4*a*b*c*z08 + 4*a*b*c*z09 - 4*a*b*c*z10 + 4*a*b*c*z11 + 2*a*b*c*z12 - 2*a*b*c*z14 + 2*a*b*c*z16 - 2*a*b*c*z18 - 2*a*b*z00 - 2*a*b*z02 + 2*a*b*z04 + 2*a*b*z06 + 4*a*b*z08 + 4*a*b*z09 - 4*a*b*z10 - 4*a*b*z11 - 2*a*b*z12 - 2*a*b*z14 + 2*a*b*z16 + 2*a*b*z18 + a*c*c*z00 - 2*a*c*c*z01 + a*c*c*z02 - a*c*c*z04 + 2*a*c*c*z05 - a*c*c*z06 - a*c*c*z12 + 2*a*c*c*z13 - a*c*c*z14 + a*c*c*z16 - 2*a*c*c*z17 + a*c*c*z18 - a*c*z00 + a*c*z02 - a*c*z04 + a*c*z06 + a*c*z12 - a*c*z14 + a*c*z16 - a*c*z18 + 2*a*z01 - 2*a*z05 - 2*a*z13 + 2*a*z17 - 2*b*c*z00 + 2*b*c*z02 + 2*b*c*z04 - 2*b*c*z06 + 4*b*c*z08 - 4*b*c*z09 - 4*b*c*z10 + 4*b*c*z11 - 2*b*c*z12 + 2*b*c*z14 + 2*b*c*z16 - 2*b*c*z18 + 2*b*z00 + 2*b*z02 + 2*b*z04 + 2*b*z06 - 4*b*z08 - 4*b*z09 - 4*b*z10 - 4*b*z11 + 2*b*z12 + 2*b*z14 + 2*b*z16 + 2*b*z18 - c*c*z00 + 2*c*c*z01 - c*c*z02 - c*c*z04 + 2*c*c*z05 - c*c*z06 + c*c*z12 - 2*c*c*z13 + c*c*z14 + c*c*z16 - 2*c*c*z17 + c*c*z18 - 2*c*z03 + 2*c*z07 + 2*c*z15 - 2*c*z19 + z00 - 2*z01 + z02 - 2*z03 + z04 - 2*z05 + z06 - 2*z07 - z12 + 2*z13 - z14 + 2*z15 - z16 + 2*z17 - z18 + 2*z19);
            double derivative_z_c = 0.125*(a*a*b*z00 - a*a*b*z02 + 2*a*a*b*z03 - a*a*b*z04 + a*a*b*z06 - 2*a*a*b*z07 - a*a*b*z12 + a*a*b*z14 - 2*a*a*b*z15 + a*a*b*z16 - a*a*b*z18 + 2*a*a*b*z19 - a*a*z00 + a*a*z02 - 2*a*a*z03 + a*a*z04 - a*a*z06 + 2*a*a*z07 - a*a*z12 + a*a*z14 - 2*a*a*z15 + a*a*z16 - a*a*z18 + 2*a*a*z19 + a*b*b*z00 - a*b*b*z02 + a*b*b*z04 - a*b*b*z06 - 2*a*b*b*z08 + 2*a*b*b*z09 - 2*a*b*b*z10 + 2*a*b*b*z11 + a*b*b*z12 - a*b*b*z14 + a*b*b*z16 - a*b*b*z18 + 2*a*b*c*z00 - 4*a*b*c*z01 + 2*a*b*c*z02 - 2*a*b*c*z04 + 4*a*b*c*z05 - 2*a*b*c*z06 - 2*a*b*c*z12 + 4*a*b*c*z13 - 2*a*b*c*z14 + 2*a*b*c*z16 - 4*a*b*c*z17 + 2*a*b*c*z18 - a*b*z00 + a*b*z02 - a*b*z04 + a*b*z06 + a*b*z12 - a*b*z14 + a*b*z16 - a*b*z18 - 2*a*c*z00 + 4*a*c*z01 - 2*a*c*z02 + 2*a*c*z04 - 4*a*c*z05 + 2*a*c*z06 - 2*a*c*z12 + 4*a*c*z13 - 2*a*c*z14 + 2*a*c*z16 - 4*a*c*z17 + 2*a*c*z18 + 2*a*z08 - 2*a*z09 + 2*a*z10 - 2*a*z11 - b*b*z00 + b*b*z02 + b*b*z04 - b*b*z06 + 2*b*b*z08 - 2*b*b*z09 - 2*b*b*z10 + 2*b*b*z11 - b*b*z12 + b*b*z14 + b*b*z16 - b*b*z18 - 2*b*c*z00 + 4*b*c*z01 - 2*b*c*z02 - 2*b*c*z04 + 4*b*c*z05 - 2*b*c*z06 + 2*b*c*z12 - 4*b*c*z13 + 2*b*c*z14 + 2*b*c*z16 - 4*b*c*z17 + 2*b*c*z18 - 2*b*z03 + 2*b*z07 + 2*b*z15 - 2*b*z19 + 2*c*z00 - 4*c*z01 + 2*c*z02 + 2*c*z04 - 4*c*z05 + 2*c*z06 + 2*c*z12 - 4*c*z13 + 2*c*z14 + 2*c*z16 - 4*c*z17 + 2*c*z18 + z00 - z02 + 2*z03 - z04 + z06 - 2*z07 - 2*z08 + 2*z09 + 2*z10 - 2*z11 + z12 - z14 + 2*z15 - z16 + z18 - 2*z19);

            // get jacobian and its inverse and determinant
            Eigen::Matrix3d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_x_c, derivative_y_a, derivative_y_b, derivative_y_c, derivative_z_a, derivative_z_b, derivative_z_c;
            Eigen::Matrix3d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int i = 0; i < 20; i++)
            {

                // get hex20 test function M
                double M = 0.;
                switch (i)
                {
                    case  0: M = 0.125*(1 - a)*(1 - b)*(1 - c)*(-a - b - c - 2); break;
                    case  1: M = 0.250*(1 - a)*(1 - b)*(1 - c*c);                break;
                    case  2: M = 0.125*(1 - a)*(1 - b)*(1 + c)*(-a - b + c - 2); break;
                    case  3: M = 0.250*(1 - a*a)*(1 - b)*(1 + c);                break;
                    case  4: M = 0.125*(1 + a)*(1 - b)*(1 + c)*(+a - b + c - 2); break;
                    case  5: M = 0.250*(1 + a)*(1 - b)*(1 - c*c);                break;
                    case  6: M = 0.125*(1 + a)*(1 - b)*(1 - c)*(+a - b - c - 2); break;
                    case  7: M = 0.250*(1 - a*a)*(1 - b)*(1 - c);                break;
                    case  8: M = 0.250*(1 - a)*(1 - b*b)*(1 - c);                break;
                    case  9: M = 0.250*(1 - a)*(1 - b*b)*(1 + c);                break;
                    case 10: M = 0.250*(1 + a)*(1 - b*b)*(1 + c);                break;
                    case 11: M = 0.250*(1 + a)*(1 - b*b)*(1 - c);                break;
                    case 12: M = 0.125*(1 - a)*(1 + b)*(1 - c)*(-a + b - c - 2); break;
                    case 13: M = 0.250*(1 - a)*(1 + b)*(1 - c*c);                break;
                    case 14: M = 0.125*(1 - a)*(1 + b)*(1 + c)*(-a + b + c - 2); break;
                    case 15: M = 0.250*(1 - a*a)*(1 + b)*(1 + c);                break;
                    case 16: M = 0.125*(1 + a)*(1 + b)*(1 + c)*(+a + b + c - 2); break;
                    case 17: M = 0.250*(1 + a)*(1 + b)*(1 - c*c);                break;
                    case 18: M = 0.125*(1 + a)*(1 + b)*(1 - c)*(+a + b - c - 2); break;
                    case 19: M = 0.250*(1 - a*a)*(1 + b)*(1 - c);                break;
                }   

                // get derivatives of hex20 test function M
                double derivative_M_a = 0.;
                double derivative_M_b = 0.;
                double derivative_M_c = 0.;
                switch (i)
                {
                    case  0: derivative_M_a = +0.125*(b - 1)*(c - 1)*(2*a + b + c + 1); derivative_M_b = +0.125*(a - 1)*(c - 1)*(a + 2*b + c + 1); derivative_M_c = +0.125*(a - 1)*(b - 1)*(a + b + 2*c + 1); break;
                    case  1: derivative_M_a = -0.250*(b - 1)*(c - 1)*(c + 1);           derivative_M_b = -0.250*(a - 1)*(c - 1)*(c + 1);           derivative_M_c = -0.500*c*(a - 1)*(b - 1);                 break;
                    case  2: derivative_M_a = -0.125*(b - 1)*(c + 1)*(2*a + b - c + 1); derivative_M_b = -0.125*(a - 1)*(c + 1)*(a + 2*b - c + 1); derivative_M_c = -0.125*(a - 1)*(b - 1)*(a + b - 2*c + 1); break;
                    case  3: derivative_M_a = +0.500*a*(b - 1)*(c + 1);                 derivative_M_b = +0.250*(a - 1)*(a + 1)*(c + 1);           derivative_M_c = +0.250*(a - 1)*(a + 1)*(b - 1);           break;
                    case  4: derivative_M_a = -0.125*(b - 1)*(c + 1)*(2*a - b + c - 1); derivative_M_b = -0.125*(a + 1)*(c + 1)*(a - 2*b + c - 1); derivative_M_c = -0.125*(a + 1)*(b - 1)*(a - b + 2*c - 1); break;
                    case  5: derivative_M_a = +0.250*(b - 1)*(c - 1)*(c + 1);           derivative_M_b = +0.250*(a + 1)*(c - 1)*(c + 1);           derivative_M_c = +0.500*c*(a + 1)*(b - 1);                 break;
                    case  6: derivative_M_a = +0.125*(b - 1)*(c - 1)*(2*a - b - c - 1); derivative_M_b = +0.125*(a + 1)*(c - 1)*(a - 2*b - c - 1); derivative_M_c = +0.125*(a + 1)*(b - 1)*(a - b - 2*c - 1); break;
                    case  7: derivative_M_a = -0.500*a*(b - 1)*(c - 1);                 derivative_M_b = -0.250*(a - 1)*(a + 1)*(c - 1);           derivative_M_c = -0.250*(a - 1)*(a + 1)*(b - 1);           break;
                    case  8: derivative_M_a = -0.250*(b - 1)*(b + 1)*(c - 1);           derivative_M_b = -0.500*b*(a - 1)*(c - 1);                 derivative_M_c = -0.250*(a - 1)*(b - 1)*(b + 1);           break;
                    case  9: derivative_M_a = +0.250*(b - 1)*(b + 1)*(c + 1);           derivative_M_b = +0.500*b*(a - 1)*(c + 1);                 derivative_M_c = +0.250*(a - 1)*(b - 1)*(b + 1);           break;
                    case 10: derivative_M_a = -0.250*(b - 1)*(b + 1)*(c + 1);           derivative_M_b = -0.500*b*(a + 1)*(c + 1);                 derivative_M_c = -0.250*(a + 1)*(b - 1)*(b + 1);           break;
                    case 11: derivative_M_a = +0.250*(b - 1)*(b + 1)*(c - 1);           derivative_M_b = +0.500*b*(a + 1)*(c - 1);                 derivative_M_c = +0.250*(a + 1)*(b - 1)*(b + 1);           break;
                    case 12: derivative_M_a = -0.125*(b + 1)*(c - 1)*(2*a - b + c + 1); derivative_M_b = -0.125*(a - 1)*(c - 1)*(a - 2*b + c + 1); derivative_M_c = -0.125*(a - 1)*(b + 1)*(a - b + 2*c + 1); break;
                    case 13: derivative_M_a = +0.250*(b + 1)*(c - 1)*(c + 1);           derivative_M_b = +0.250*(a - 1)*(c - 1)*(c + 1);           derivative_M_c = +0.500*c*(a - 1)*(b + 1);                 break;
                    case 14: derivative_M_a = +0.125*(b + 1)*(c + 1)*(2*a - b - c + 1); derivative_M_b = +0.125*(a - 1)*(c + 1)*(a - 2*b - c + 1); derivative_M_c = +0.125*(a - 1)*(b + 1)*(a - b - 2*c + 1); break;
                    case 15: derivative_M_a = -0.500*a*(b + 1)*(c + 1);                 derivative_M_b = -0.250*(a - 1)*(a + 1)*(c + 1);           derivative_M_c = -0.250*(a - 1)*(a + 1)*(b + 1);           break;
                    case 16: derivative_M_a = +0.125*(b + 1)*(c + 1)*(2*a + b + c - 1); derivative_M_b = +0.125*(a + 1)*(c + 1)*(a + 2*b + c - 1); derivative_M_c = +0.125*(a + 1)*(b + 1)*(a + b + 2*c - 1); break;
                    case 17: derivative_M_a = -0.250*(b + 1)*(c - 1)*(c + 1);           derivative_M_b = -0.250*(a + 1)*(c - 1)*(c + 1);           derivative_M_c = -0.500*c*(a + 1)*(b + 1);                 break;
                    case 18: derivative_M_a = -0.125*(b + 1)*(c - 1)*(2*a + b - c - 1); derivative_M_b = -0.125*(a + 1)*(c - 1)*(a + 2*b - c - 1); derivative_M_c = -0.125*(a + 1)*(b + 1)*(a + b - 2*c - 1); break;
                    case 19: derivative_M_a = +0.500*a*(b + 1)*(c - 1);                 derivative_M_b = +0.250*(a - 1)*(a + 1)*(c - 1);           derivative_M_c = +0.250*(a - 1)*(a + 1)*(b + 1);           break;
                }

                // get vectors with derivatives of test functions wrt a, b, and c
                Eigen::RowVector3d derivative_M_abc_vec;
                derivative_M_abc_vec << derivative_M_a, derivative_M_b, derivative_M_c;

                // get vectors with derivatives of test functions wrt x, y, and z
                Eigen::RowVector3d derivative_M_xyz_vec;
                derivative_M_xyz_vec = derivative_M_abc_vec*jacobian_inverse_mat;
                double derivative_M_x = derivative_M_xyz_vec.coeffRef(0);
                double derivative_M_y = derivative_M_xyz_vec.coeffRef(1);
                double derivative_M_z = derivative_M_xyz_vec.coeffRef(2);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_M_x_part_mli_vec.push_back(derivative_M_x);
                derivative_M_y_part_mli_vec.push_back(derivative_M_y);
                derivative_M_z_part_mli_vec.push_back(derivative_M_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_M_x_part_ml_vec.push_back(derivative_M_x_part_mli_vec);
            derivative_M_y_part_ml_vec.push_back(derivative_M_y_part_mli_vec);
            derivative_M_z_part_ml_vec.push_back(derivative_M_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        M_vec.push_back(M_part_ml_vec);
        derivative_M_x_vec.push_back(derivative_M_x_part_ml_vec);
        derivative_M_y_vec.push_back(derivative_M_y_part_ml_vec);
        derivative_M_z_vec.push_back(derivative_M_z_part_ml_vec);

    }

    // iterate for each hex8 element
    for (int m = 0; m < gh8s.num_element; m++)
    {

        // initialize
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
        for (int l = 0; l < 27; l++)
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
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);
            derivative_N_y_part_ml_vec.push_back(derivative_N_y_part_mli_vec);
            derivative_N_z_part_ml_vec.push_back(derivative_N_z_part_mli_vec);

        }

        // store in vectors
        N_vec.push_back(N_part_ml_vec);
        derivative_N_x_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_N_y_vec.push_back(derivative_N_y_part_ml_vec);
        derivative_N_z_vec.push_back(derivative_N_z_part_ml_vec);
        
    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<double> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i];
        }

    integral_part_i_vec.push_back(integral_value);
    }
    integral_Mi_hex20_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Ni_hex8_derivative_Mj_hex20_x()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * derivative_M_x_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_hex8_derivative_Mj_hex20_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Ni_hex8_derivative_Mj_hex20_y()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * derivative_M_y_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_hex8_derivative_Mj_hex20_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Ni_hex8_derivative_Mj_hex20_z()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 8; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * N_vec[m][l][i] * derivative_M_z_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_hex8_derivative_Mj_hex20_z_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_derivative_Nj_hex8_x()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 8; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * derivative_N_x_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_derivative_Nj_hex8_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_derivative_Nj_hex8_y()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 8; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * derivative_N_y_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_derivative_Nj_hex8_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_derivative_Nj_hex8_z()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 8; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * derivative_N_z_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_derivative_Nj_hex8_z_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_div_Mi_hex20_dot_div_Mj_hex20()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * (derivative_M_x_vec[m][l][i]*derivative_M_x_vec[m][l][j] + derivative_M_y_vec[m][l][i]*derivative_M_y_vec[m][l][j] + derivative_M_z_vec[m][l][i]*derivative_M_z_vec[m][l][j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Mi_hex20_dot_div_Mj_hex20_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_Mj_hex20()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<double>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<double> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * M_vec[m][l][j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_Mj_hex20_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<std::vector<double>>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<std::vector<double>> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){  
    std::vector<double> integral_part_ijk_vec;
    for (int k = 0; k < 20; k++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * M_vec[m][l][k] * derivative_M_x_vec[m][l][j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<std::vector<double>>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<std::vector<double>> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){  
    std::vector<double> integral_part_ijk_vec;
    for (int k = 0; k < 20; k++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * M_vec[m][l][k] * derivative_M_y_vec[m][l][j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y_vec.push_back(integral_part_i_vec);

    }

}

void IntegralHex8Hex20Class::evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z()
{
    
    // iterate for each element
    for (int m = 0; m < gh20s.num_element; m++){  
    
    // iterate for each test function combination
    std::vector<std::vector<std::vector<double>>> integral_part_i_vec;
    for (int i = 0; i < 20; i++){  
    std::vector<std::vector<double>> integral_part_ij_vec;
    for (int j = 0; j < 20; j++){  
    std::vector<double> integral_part_ijk_vec;
    for (int k = 0; k < 20; k++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int l = 0; l < 27; l++) 
        {
            integral_value += jacobian_determinant_vec[m][l] * M_vec[m][l][i] * M_vec[m][l][k] * derivative_M_z_vec[m][l][j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z_vec.push_back(integral_part_i_vec);

    }

}

#endif
