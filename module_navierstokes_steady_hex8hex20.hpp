#ifndef MODULE_NAVIERSTOKES_STEADY_HEX8HEX20
#define MODULE_NAVIERSTOKES_STEADY_HEX8HEX20
#include <vector>
#include "Eigen/Eigen"
#include "integral_hex8hex20.hpp"
#include "grid_hex8.hpp"
#include "grid_hex20.hpp"
#include "scalar_hex20.hpp"

class NavierStokesSteadyHex8Hex20Class
{

    public:

    // variables
    GridHex20Struct gh20s;
    GridHex8Struct gh8s;
    BoundaryHex20Struct u_bh20s;
    BoundaryHex20Struct v_bh20s;
    BoundaryHex20Struct w_bh20s;
    BoundaryHex8Struct p_bh8s;
    IntegralHex8Hex20Class ih8h20c;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
        ScalarHex20Class &density_sh20c, ScalarHex20Class &viscosity_sh20c, double gravity_x, double gravity_y, double gravity_z,
        int start_id_u, int start_id_v, int start_id_w, int start_id_p
    );

    // constructor
    NavierStokesSteadyHex8Hex20Class(
        GridHex20Struct &gh20s_in, GridHex8Struct &gh8s_in,
        BoundaryHex20Struct &u_bh20s_in, BoundaryHex20Struct &v_bh20s_in, BoundaryHex20Struct &w_bh20s_in, BoundaryHex8Struct &p_bh8s_in,
        IntegralHex8Hex20Class &ih8h20c_in
    )
    {
        
        // variables
        gh20s = gh20s_in;
        gh8s = gh8s_in;
        u_bh20s = u_bh20s_in;
        v_bh20s = v_bh20s_in;
        w_bh20s = w_bh20s_in;
        p_bh8s = p_bh8s_in;
        ih8h20c = ih8h20c_in;

        // integrals
        ih8h20c.evaluate_test_functions_derivatives();
        ih8h20c.evaluate_integral_Mi_hex20();
        ih8h20c.evaluate_integral_Ni_hex8_derivative_Mj_hex20_x();
        ih8h20c.evaluate_integral_Ni_hex8_derivative_Mj_hex20_y();
        ih8h20c.evaluate_integral_Ni_hex8_derivative_Mj_hex20_z();
        ih8h20c.evaluate_integral_Mi_hex20_derivative_Nj_hex8_x();
        ih8h20c.evaluate_integral_Mi_hex20_derivative_Nj_hex8_y();
        ih8h20c.evaluate_integral_Mi_hex20_derivative_Nj_hex8_z();
        ih8h20c.evaluate_integral_div_Mi_hex20_dot_div_Mj_hex20();
        ih8h20c.evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x();
        ih8h20c.evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y();
        ih8h20c.evaluate_integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z();

    }

};

void NavierStokesSteadyHex8Hex20Class::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_last_iteration_vec,
    ScalarHex20Class &density_sh20c, ScalarHex20Class &viscosity_sh20c, double gravity_x, double gravity_y, double gravity_z,
    int start_id_u, int start_id_v, int start_id_w, int start_id_p
)
{

    // fill up matrix with NSE
    for (int m = 0; m < gh20s.num_element; m++)
    {

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get id of points around hex8 element
        int n00_hex8 = gh8s.element_p00_id_vec[m]; int n04_hex8 = gh8s.element_p04_id_vec[m];
        int n01_hex8 = gh8s.element_p01_id_vec[m]; int n05_hex8 = gh8s.element_p05_id_vec[m];
        int n02_hex8 = gh8s.element_p02_id_vec[m]; int n06_hex8 = gh8s.element_p06_id_vec[m];
        int n03_hex8 = gh8s.element_p03_id_vec[m]; int n07_hex8 = gh8s.element_p07_id_vec[m];
        int n_hex8_arr[8] = {n00_hex8, n01_hex8, n02_hex8, n03_hex8, n04_hex8, n05_hex8, n06_hex8, n07_hex8};

        // get u velocity values at points around element
        double u00_hex20 = x_last_iteration_vec[start_id_u + n00_hex20]; double u10_hex20 = x_last_iteration_vec[start_id_u + n10_hex20];
        double u01_hex20 = x_last_iteration_vec[start_id_u + n01_hex20]; double u11_hex20 = x_last_iteration_vec[start_id_u + n11_hex20];
        double u02_hex20 = x_last_iteration_vec[start_id_u + n02_hex20]; double u12_hex20 = x_last_iteration_vec[start_id_u + n12_hex20];
        double u03_hex20 = x_last_iteration_vec[start_id_u + n03_hex20]; double u13_hex20 = x_last_iteration_vec[start_id_u + n13_hex20];
        double u04_hex20 = x_last_iteration_vec[start_id_u + n04_hex20]; double u14_hex20 = x_last_iteration_vec[start_id_u + n14_hex20];
        double u05_hex20 = x_last_iteration_vec[start_id_u + n05_hex20]; double u15_hex20 = x_last_iteration_vec[start_id_u + n15_hex20];
        double u06_hex20 = x_last_iteration_vec[start_id_u + n06_hex20]; double u16_hex20 = x_last_iteration_vec[start_id_u + n16_hex20];
        double u07_hex20 = x_last_iteration_vec[start_id_u + n07_hex20]; double u17_hex20 = x_last_iteration_vec[start_id_u + n17_hex20];
        double u08_hex20 = x_last_iteration_vec[start_id_u + n08_hex20]; double u18_hex20 = x_last_iteration_vec[start_id_u + n18_hex20];
        double u09_hex20 = x_last_iteration_vec[start_id_u + n09_hex20]; double u19_hex20 = x_last_iteration_vec[start_id_u + n19_hex20];
        double u_hex20_arr[20] = {
            u00_hex20, u01_hex20, u02_hex20, u03_hex20, u04_hex20, u05_hex20, u06_hex20, u07_hex20, u08_hex20, u09_hex20,
            u10_hex20, u11_hex20, u12_hex20, u13_hex20, u14_hex20, u15_hex20, u16_hex20, u17_hex20, u18_hex20, u19_hex20
        };

        // get v velocity values at points around element
        double v00_hex20 = x_last_iteration_vec[start_id_v + n00_hex20]; double v10_hex20 = x_last_iteration_vec[start_id_v + n10_hex20];
        double v01_hex20 = x_last_iteration_vec[start_id_v + n01_hex20]; double v11_hex20 = x_last_iteration_vec[start_id_v + n11_hex20];
        double v02_hex20 = x_last_iteration_vec[start_id_v + n02_hex20]; double v12_hex20 = x_last_iteration_vec[start_id_v + n12_hex20];
        double v03_hex20 = x_last_iteration_vec[start_id_v + n03_hex20]; double v13_hex20 = x_last_iteration_vec[start_id_v + n13_hex20];
        double v04_hex20 = x_last_iteration_vec[start_id_v + n04_hex20]; double v14_hex20 = x_last_iteration_vec[start_id_v + n14_hex20];
        double v05_hex20 = x_last_iteration_vec[start_id_v + n05_hex20]; double v15_hex20 = x_last_iteration_vec[start_id_v + n15_hex20];
        double v06_hex20 = x_last_iteration_vec[start_id_v + n06_hex20]; double v16_hex20 = x_last_iteration_vec[start_id_v + n16_hex20];
        double v07_hex20 = x_last_iteration_vec[start_id_v + n07_hex20]; double v17_hex20 = x_last_iteration_vec[start_id_v + n17_hex20];
        double v08_hex20 = x_last_iteration_vec[start_id_v + n08_hex20]; double v18_hex20 = x_last_iteration_vec[start_id_v + n18_hex20];
        double v09_hex20 = x_last_iteration_vec[start_id_v + n09_hex20]; double v19_hex20 = x_last_iteration_vec[start_id_v + n19_hex20];
        double v_hex20_arr[20] = {
            v00_hex20, v01_hex20, v02_hex20, v03_hex20, v04_hex20, v05_hex20, v06_hex20, v07_hex20, v08_hex20, v09_hex20,
            v10_hex20, v11_hex20, v12_hex20, v13_hex20, v14_hex20, v15_hex20, v16_hex20, v17_hex20, v18_hex20, v19_hex20
        };

        // get w velocity values at points around element
        double w00_hex20 = x_last_iteration_vec[start_id_w + n00_hex20]; double w10_hex20 = x_last_iteration_vec[start_id_w + n10_hex20];
        double w01_hex20 = x_last_iteration_vec[start_id_w + n01_hex20]; double w11_hex20 = x_last_iteration_vec[start_id_w + n11_hex20];
        double w02_hex20 = x_last_iteration_vec[start_id_w + n02_hex20]; double w12_hex20 = x_last_iteration_vec[start_id_w + n12_hex20];
        double w03_hex20 = x_last_iteration_vec[start_id_w + n03_hex20]; double w13_hex20 = x_last_iteration_vec[start_id_w + n13_hex20];
        double w04_hex20 = x_last_iteration_vec[start_id_w + n04_hex20]; double w14_hex20 = x_last_iteration_vec[start_id_w + n14_hex20];
        double w05_hex20 = x_last_iteration_vec[start_id_w + n05_hex20]; double w15_hex20 = x_last_iteration_vec[start_id_w + n15_hex20];
        double w06_hex20 = x_last_iteration_vec[start_id_w + n06_hex20]; double w16_hex20 = x_last_iteration_vec[start_id_w + n16_hex20];
        double w07_hex20 = x_last_iteration_vec[start_id_w + n07_hex20]; double w17_hex20 = x_last_iteration_vec[start_id_w + n17_hex20];
        double w08_hex20 = x_last_iteration_vec[start_id_w + n08_hex20]; double w18_hex20 = x_last_iteration_vec[start_id_w + n18_hex20];
        double w09_hex20 = x_last_iteration_vec[start_id_w + n09_hex20]; double w19_hex20 = x_last_iteration_vec[start_id_w + n19_hex20];
        double w_hex20_arr[20] = {
            w00_hex20, w01_hex20, w02_hex20, w03_hex20, w04_hex20, w05_hex20, w06_hex20, w07_hex20, w08_hex20, w09_hex20,
            w10_hex20, w11_hex20, w12_hex20, w13_hex20, w14_hex20, w15_hex20, w16_hex20, w17_hex20, w18_hex20, w19_hex20
        };

        // get density at points around element
        double density00_hex20 = density_sh20c.scalar_vec[n00_hex20]; double density10_hex20 = density_sh20c.scalar_vec[n10_hex20];
        double density01_hex20 = density_sh20c.scalar_vec[n01_hex20]; double density11_hex20 = density_sh20c.scalar_vec[n11_hex20];
        double density02_hex20 = density_sh20c.scalar_vec[n02_hex20]; double density12_hex20 = density_sh20c.scalar_vec[n12_hex20];
        double density03_hex20 = density_sh20c.scalar_vec[n03_hex20]; double density13_hex20 = density_sh20c.scalar_vec[n13_hex20];
        double density04_hex20 = density_sh20c.scalar_vec[n04_hex20]; double density14_hex20 = density_sh20c.scalar_vec[n14_hex20];
        double density05_hex20 = density_sh20c.scalar_vec[n05_hex20]; double density15_hex20 = density_sh20c.scalar_vec[n15_hex20];
        double density06_hex20 = density_sh20c.scalar_vec[n06_hex20]; double density16_hex20 = density_sh20c.scalar_vec[n16_hex20];
        double density07_hex20 = density_sh20c.scalar_vec[n07_hex20]; double density17_hex20 = density_sh20c.scalar_vec[n17_hex20];
        double density08_hex20 = density_sh20c.scalar_vec[n08_hex20]; double density18_hex20 = density_sh20c.scalar_vec[n18_hex20];
        double density09_hex20 = density_sh20c.scalar_vec[n09_hex20]; double density19_hex20 = density_sh20c.scalar_vec[n19_hex20];
        double density_hex20_arr[20] = {
            density00_hex20, density01_hex20, density02_hex20, density03_hex20, density04_hex20, density05_hex20, density06_hex20, density07_hex20, density08_hex20, density09_hex20,
            density10_hex20, density11_hex20, density12_hex20, density13_hex20, density14_hex20, density15_hex20, density16_hex20, density17_hex20, density18_hex20, density19_hex20
        };

        // get viscosity at points around element
        double viscosity00_hex20 = viscosity_sh20c.scalar_vec[n00_hex20]; double viscosity10_hex20 = viscosity_sh20c.scalar_vec[n10_hex20];
        double viscosity01_hex20 = viscosity_sh20c.scalar_vec[n01_hex20]; double viscosity11_hex20 = viscosity_sh20c.scalar_vec[n11_hex20];
        double viscosity02_hex20 = viscosity_sh20c.scalar_vec[n02_hex20]; double viscosity12_hex20 = viscosity_sh20c.scalar_vec[n12_hex20];
        double viscosity03_hex20 = viscosity_sh20c.scalar_vec[n03_hex20]; double viscosity13_hex20 = viscosity_sh20c.scalar_vec[n13_hex20];
        double viscosity04_hex20 = viscosity_sh20c.scalar_vec[n04_hex20]; double viscosity14_hex20 = viscosity_sh20c.scalar_vec[n14_hex20];
        double viscosity05_hex20 = viscosity_sh20c.scalar_vec[n05_hex20]; double viscosity15_hex20 = viscosity_sh20c.scalar_vec[n15_hex20];
        double viscosity06_hex20 = viscosity_sh20c.scalar_vec[n06_hex20]; double viscosity16_hex20 = viscosity_sh20c.scalar_vec[n16_hex20];
        double viscosity07_hex20 = viscosity_sh20c.scalar_vec[n07_hex20]; double viscosity17_hex20 = viscosity_sh20c.scalar_vec[n17_hex20];
        double viscosity08_hex20 = viscosity_sh20c.scalar_vec[n08_hex20]; double viscosity18_hex20 = viscosity_sh20c.scalar_vec[n18_hex20];
        double viscosity09_hex20 = viscosity_sh20c.scalar_vec[n09_hex20]; double viscosity19_hex20 = viscosity_sh20c.scalar_vec[n19_hex20];
        double viscosity_hex20_arr[20] = {
            viscosity00_hex20, viscosity01_hex20, viscosity02_hex20, viscosity03_hex20, viscosity04_hex20, viscosity05_hex20, viscosity06_hex20, viscosity07_hex20, viscosity08_hex20, viscosity09_hex20,
            viscosity10_hex20, viscosity11_hex20, viscosity12_hex20, viscosity13_hex20, viscosity14_hex20, viscosity15_hex20, viscosity16_hex20, viscosity17_hex20, viscosity18_hex20, viscosity19_hex20
        };

        // NSE x-component

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {

            // calculate velocity integrals
            double integral_Mi_hex20_u_derivative_Mj_hex20_x = 0.;
            double integral_Mi_hex20_v_derivative_Mj_hex20_y = 0.;
            double integral_Mi_hex20_w_derivative_Mj_hex20_z = 0.;
            for (int k = 0; k < 20; k++)
            {
                integral_Mi_hex20_u_derivative_Mj_hex20_x += u_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x_vec[m][i][j][k];
                integral_Mi_hex20_v_derivative_Mj_hex20_y += v_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y_vec[m][i][j][k];
                integral_Mi_hex20_w_derivative_Mj_hex20_z += w_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z_vec[m][i][j][k];
            }

            // fill up coefficient matrix
            a_mat.coeffRef(start_id_u + n_hex20_arr[i], start_id_u + n_hex20_arr[j]) += (
                  viscosity_hex20_arr[i]*ih8h20c.integral_div_Mi_hex20_dot_div_Mj_hex20_vec[m][i][j]
                + density_hex20_arr[i]*integral_Mi_hex20_u_derivative_Mj_hex20_x
                + density_hex20_arr[i]*integral_Mi_hex20_v_derivative_Mj_hex20_y
                + density_hex20_arr[i]*integral_Mi_hex20_w_derivative_Mj_hex20_z
            );

        }}

        // coefficients of p
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_u + n_hex20_arr[i], start_id_p + n_hex8_arr[j]) += ih8h20c.integral_Mi_hex20_derivative_Nj_hex8_x_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 20; i++) {
            b_vec.coeffRef(start_id_u + n_hex20_arr[i]) += density_hex20_arr[i]*gravity_x*ih8h20c.integral_Mi_hex20_vec[m][i];
        }

        // NSE y-component

        // coefficients of v
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {

            // calculate velocity integrals
            double integral_Mi_hex20_u_derivative_Mj_hex20_x = 0.;
            double integral_Mi_hex20_v_derivative_Mj_hex20_y = 0.;
            double integral_Mi_hex20_w_derivative_Mj_hex20_z = 0.;
            for (int k = 0; k < 20; k++)
            {
                integral_Mi_hex20_u_derivative_Mj_hex20_x += u_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x_vec[m][i][j][k];
                integral_Mi_hex20_v_derivative_Mj_hex20_y += v_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y_vec[m][i][j][k];
                integral_Mi_hex20_w_derivative_Mj_hex20_z += w_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z_vec[m][i][j][k];
            }

            // fill up coefficient matrix
            a_mat.coeffRef(start_id_v + n_hex20_arr[i], start_id_v + n_hex20_arr[j]) += (
                  viscosity_hex20_arr[i]*ih8h20c.integral_div_Mi_hex20_dot_div_Mj_hex20_vec[m][i][j]
                + density_hex20_arr[i]*integral_Mi_hex20_u_derivative_Mj_hex20_x
                + density_hex20_arr[i]*integral_Mi_hex20_v_derivative_Mj_hex20_y
                + density_hex20_arr[i]*integral_Mi_hex20_w_derivative_Mj_hex20_z
            );

        }}

        // coefficients of p
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_v + n_hex20_arr[i], start_id_p + n_hex8_arr[j]) += ih8h20c.integral_Mi_hex20_derivative_Nj_hex8_y_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 20; i++) {
            b_vec.coeffRef(start_id_v + n_hex20_arr[i]) += density_hex20_arr[i]*gravity_y*ih8h20c.integral_Mi_hex20_vec[m][i];
        }

        // NSE z-component

        // coefficients of w
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {

            // calculate velocity integrals
            double integral_Mi_hex20_u_derivative_Mj_hex20_x = 0.;
            double integral_Mi_hex20_v_derivative_Mj_hex20_y = 0.;
            double integral_Mi_hex20_w_derivative_Mj_hex20_z = 0.;
            for (int k = 0; k < 20; k++)
            {
                integral_Mi_hex20_u_derivative_Mj_hex20_x += u_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_x_vec[m][i][j][k];
                integral_Mi_hex20_v_derivative_Mj_hex20_y += v_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_y_vec[m][i][j][k];
                integral_Mi_hex20_w_derivative_Mj_hex20_z += w_hex20_arr[k] * ih8h20c.integral_Mi_hex20_Mk_hex20_derivative_Mj_hex20_z_vec[m][i][j][k];
            }

            // fill up coefficient matrix
            a_mat.coeffRef(start_id_w + n_hex20_arr[i], start_id_w + n_hex20_arr[j]) += (
                  viscosity_hex20_arr[i]*ih8h20c.integral_div_Mi_hex20_dot_div_Mj_hex20_vec[m][i][j]
                + density_hex20_arr[i]*integral_Mi_hex20_u_derivative_Mj_hex20_x
                + density_hex20_arr[i]*integral_Mi_hex20_v_derivative_Mj_hex20_y
                + density_hex20_arr[i]*integral_Mi_hex20_w_derivative_Mj_hex20_z
            );

        }}

        // coefficients of p
        for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 8; j++) {
            a_mat.coeffRef(start_id_w + n_hex20_arr[i], start_id_p + n_hex8_arr[j]) += ih8h20c.integral_Mi_hex20_derivative_Nj_hex8_z_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 20; i++) {
            b_vec.coeffRef(start_id_w + n_hex20_arr[i]) += density_hex20_arr[i]*gravity_z*ih8h20c.integral_Mi_hex20_vec[m][i];
        }

    }

    // fill up matrix with continuity equation
    for (int m = 0; m < gh8s.num_element; m++)
    {

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get id of points around hex8 element
        int n00_hex8 = gh8s.element_p00_id_vec[m]; int n04_hex8 = gh8s.element_p04_id_vec[m];
        int n01_hex8 = gh8s.element_p01_id_vec[m]; int n05_hex8 = gh8s.element_p05_id_vec[m];
        int n02_hex8 = gh8s.element_p02_id_vec[m]; int n06_hex8 = gh8s.element_p06_id_vec[m];
        int n03_hex8 = gh8s.element_p03_id_vec[m]; int n07_hex8 = gh8s.element_p07_id_vec[m];
        int n_hex8_arr[8] = {n00_hex8, n01_hex8, n02_hex8, n03_hex8, n04_hex8, n05_hex8, n06_hex8, n07_hex8};

        // calculate a_mat coefficients
        // row i - one linearized equation for each test function
        // column j - coefficient of each variable [u0 ... uN v0 ... vN p0 ... pN]

        // coefficients of u
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 20; j++) {
            a_mat.coeffRef(start_id_p + n_hex8_arr[i], start_id_u + n_hex20_arr[j]) += ih8h20c.integral_Ni_hex8_derivative_Mj_hex20_x_vec[m][i][j];
        }}

        // coefficients of v
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 20; j++) {
            a_mat.coeffRef(start_id_p + n_hex8_arr[i], start_id_v + n_hex20_arr[j]) += ih8h20c.integral_Ni_hex8_derivative_Mj_hex20_y_vec[m][i][j];
        }}

        // coefficients of w
        for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 20; j++) {
            a_mat.coeffRef(start_id_p + n_hex8_arr[i], start_id_w + n_hex20_arr[j]) += ih8h20c.integral_Ni_hex8_derivative_Mj_hex20_z_vec[m][i][j];
        }}

        // calculate b_vec coefficients
        for (int i = 0; i < 8; i++) {
            b_vec.coeffRef(start_id_p + n_hex8_arr[i]) += 0.;
        }

    }

//    // iterate for each flux boundary element (x-component of NSE)
//    for (int k = 0; k < u_bh20s.num_element_flux; k++)
//    {
//
//        // get id of element
//        int m = u_bh20s.element_flux_id_vec[k];
//
//        // get id of points around hex20 element
//        int n0_hex20 = gh20s.element_p0_id_vec[m];
//        int n1_hex20 = gh20s.element_p1_id_vec[m];
//        int n2_hex20 = gh20s.element_p2_id_vec[m];
//        int n3_hex20 = gh20s.element_p3_id_vec[m];
//        int n4_hex20 = gh20s.element_p4_id_vec[m];
//        int n5_hex20 = gh20s.element_p5_id_vec[m];
//        int n6_hex20 = gh20s.element_p6_id_vec[m];
//        int n7_hex20 = gh20s.element_p7_id_vec[m];
//        int n_hex20_arr[8] = {n0_hex20, n1_hex20, n2_hex20, n3_hex20, n4_hex20, n5_hex20, n6_hex20, n7_hex20};
//
//        // get x coordinates of points around hex20 element
//        double x0_hex20 = gh20s.point_pos_x_vec[n0_hex20];
//        double x1_hex20 = gh20s.point_pos_x_vec[n1_hex20];
//        double x2_hex20 = gh20s.point_pos_x_vec[n2_hex20];
//        double x3_hex20 = gh20s.point_pos_x_vec[n3_hex20];
//        double x4_hex20 = gh20s.point_pos_x_vec[n4_hex20];
//        double x5_hex20 = gh20s.point_pos_x_vec[n5_hex20];
//        double x6_hex20 = gh20s.point_pos_x_vec[n6_hex20];
//        double x7_hex20 = gh20s.point_pos_x_vec[n7_hex20];
//        double x_hex20_arr[8] = {x0_hex20, x1_hex20, x2_hex20, x3_hex20, x4_hex20, x5_hex20, x6_hex20, x7_hex20};
//
//        // get y coordinates of points around hex20 element
//        double y0_hex20 = gh20s.point_pos_y_vec[n0_hex20];
//        double y1_hex20 = gh20s.point_pos_y_vec[n1_hex20];
//        double y2_hex20 = gh20s.point_pos_y_vec[n2_hex20];
//        double y3_hex20 = gh20s.point_pos_y_vec[n3_hex20];
//        double y4_hex20 = gh20s.point_pos_y_vec[n4_hex20];
//        double y5_hex20 = gh20s.point_pos_y_vec[n5_hex20];
//        double y6_hex20 = gh20s.point_pos_y_vec[n6_hex20];
//        double y7_hex20 = gh20s.point_pos_y_vec[n7_hex20];
//        double y_hex20_arr[8] = {y0_hex20, y1_hex20, y2_hex20, y3_hex20, y4_hex20, y5_hex20, y6_hex20, y7_hex20};
//
//        // get points where the boundary is applied
//        int a = u_bh20s.element_flux_pa_loc_vec[k];  // 0 to 7
//        int b = u_bh20s.element_flux_pb_loc_vec[k];  // 0 to 7
//        int c = u_bh20s.element_flux_pb_loc_vec[k];  // 0 to 7
//
//        // identify boundary type
//        int config_id = u_bh20s.element_flux_config_id_vec[k];
//        BoundaryConfigHex20Struct bch20s = u_bh20s.boundary_config_vec[config_id];
//
//        // apply boundary condition
//        if (bch20s.boundary_type_str == "neumann")
//        {
//
//            // calculate distance between a-b-c
//            // APPROXIMATION ONLY. FOR REVISION.
//            double xa = x_hex20_arr[a];
//            double xc = x_hex20_arr[c];
//            double ya = y_hex20_arr[a];
//            double yc = y_hex20_arr[c];
//            double dist_pa_pb_pc = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));
//
//            // add to b_vec
//            // NOT YET VERIFIED. PLACEHOLDER ONLY.
//            b_vec.coeffRef(start_id_u + n_hex20_arr[a]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//            b_vec.coeffRef(start_id_u + n_hex20_arr[b]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//            b_vec.coeffRef(start_id_u + n_hex20_arr[c]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//
//        }
//
//    }
//
//    // iterate for each flux boundary element (y-component of NSE)
//    for (int k = 0; k < v_bh20s.num_element_flux; k++)
//    {
//
//        // get id of element
//        int m = v_bh20s.element_flux_id_vec[k];
//
//        // get id of points around hex20 element
//        int n0_hex20 = gh20s.element_p0_id_vec[m];
//        int n1_hex20 = gh20s.element_p1_id_vec[m];
//        int n2_hex20 = gh20s.element_p2_id_vec[m];
//        int n3_hex20 = gh20s.element_p3_id_vec[m];
//        int n4_hex20 = gh20s.element_p4_id_vec[m];
//        int n5_hex20 = gh20s.element_p5_id_vec[m];
//        int n6_hex20 = gh20s.element_p6_id_vec[m];
//        int n7_hex20 = gh20s.element_p7_id_vec[m];
//        int n_hex20_arr[8] = {n0_hex20, n1_hex20, n2_hex20, n3_hex20, n4_hex20, n5_hex20, n6_hex20, n7_hex20};
//
//        // get x coordinates of points around hex20 element
//        double x0_hex20 = gh20s.point_pos_x_vec[n0_hex20];
//        double x1_hex20 = gh20s.point_pos_x_vec[n1_hex20];
//        double x2_hex20 = gh20s.point_pos_x_vec[n2_hex20];
//        double x3_hex20 = gh20s.point_pos_x_vec[n3_hex20];
//        double x4_hex20 = gh20s.point_pos_x_vec[n4_hex20];
//        double x5_hex20 = gh20s.point_pos_x_vec[n5_hex20];
//        double x6_hex20 = gh20s.point_pos_x_vec[n6_hex20];
//        double x7_hex20 = gh20s.point_pos_x_vec[n7_hex20];
//        double x_hex20_arr[8] = {x0_hex20, x1_hex20, x2_hex20, x3_hex20, x4_hex20, x5_hex20, x6_hex20, x7_hex20};
//
//        // get y coordinates of points around hex20 element
//        double y0_hex20 = gh20s.point_pos_y_vec[n0_hex20];
//        double y1_hex20 = gh20s.point_pos_y_vec[n1_hex20];
//        double y2_hex20 = gh20s.point_pos_y_vec[n2_hex20];
//        double y3_hex20 = gh20s.point_pos_y_vec[n3_hex20];
//        double y4_hex20 = gh20s.point_pos_y_vec[n4_hex20];
//        double y5_hex20 = gh20s.point_pos_y_vec[n5_hex20];
//        double y6_hex20 = gh20s.point_pos_y_vec[n6_hex20];
//        double y7_hex20 = gh20s.point_pos_y_vec[n7_hex20];
//        double y_hex20_arr[8] = {y0_hex20, y1_hex20, y2_hex20, y3_hex20, y4_hex20, y5_hex20, y6_hex20, y7_hex20};
//
//        // get points where the boundary is applied
//        int a = v_bh20s.element_flux_pa_loc_vec[k];  // 0 to 7
//        int b = v_bh20s.element_flux_pb_loc_vec[k];  // 0 to 7
//        int c = v_bh20s.element_flux_pb_loc_vec[k];  // 0 to 7
//
//        // identify boundary type
//        int config_id = v_bh20s.element_flux_config_id_vec[k];
//        BoundaryConfigHex20Struct bch20s = v_bh20s.boundary_config_vec[config_id];
//
//        // apply boundary condition
//        if (bch20s.boundary_type_str == "neumann")
//        {
//
//            // calculate distance between a-b-c
//            // APPROXIMATION ONLY. FOR REVISION.
//            double xa = x_hex20_arr[a];
//            double xc = x_hex20_arr[c];
//            double ya = y_hex20_arr[a];
//            double yc = y_hex20_arr[c];
//            double dist_pa_pb_pc = sqrt((xa - xc)*(xa - xc) + (ya - yc)*(ya - yc));
//
//            // add to b_vec
//            // NOT YET VERIFIED. PLACEHOLDER ONLY.
//            b_vec.coeffRef(start_id_v + n_hex20_arr[a]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//            b_vec.coeffRef(start_id_v + n_hex20_arr[b]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//            b_vec.coeffRef(start_id_v + n_hex20_arr[c]) += 0.5 * bch20s.boundary_parameter_vec[0] * dist_pa_pb_pc;
//
//        }
//
//    }
//
//    // ADD Z COMPONENT
//
//    // iterate for each flux boundary element (continuity)
//    for (int k = 0; k < p_bh8s.num_element_flux; k++)
//    {
//
//        // get id of element
//        int m = p_bh8s.element_flux_id_vec[k];
//
//        // get id of points around hex8 element
//        int n0_hex8 = gh8s.element_p0_id_vec[m];
//        int n1_hex8 = gh8s.element_p1_id_vec[m];
//        int n2_hex8 = gh8s.element_p2_id_vec[m];
//        int n3_hex8 = gh8s.element_p3_id_vec[m];
//        int n_hex8_arr[4] = {n0_hex8, n1_hex8, n2_hex8, n3_hex8};
//
//        // get x coordinates of points around element
//        double x0_hex8 = gh8s.point_pos_x_vec[n0_hex8];
//        double x1_hex8 = gh8s.point_pos_x_vec[n1_hex8];
//        double x2_hex8 = gh8s.point_pos_x_vec[n2_hex8];
//        double x3_hex8 = gh8s.point_pos_x_vec[n3_hex8];
//        double x_hex8_arr[4] = {x0_hex8, x1_hex8, x2_hex8, x3_hex8};
//
//        // get y coordinates of points around element
//        double y0_hex8 = gh8s.point_pos_y_vec[n0_hex8];
//        double y1_hex8 = gh8s.point_pos_y_vec[n1_hex8];
//        double y2_hex8 = gh8s.point_pos_y_vec[n2_hex8];
//        double y3_hex8 = gh8s.point_pos_y_vec[n3_hex8];
//        double y_hex8_arr[4] = {y0_hex8, y1_hex8, y2_hex8, y3_hex8};
//
//        // get points where the boundary is applied
//        int a = p_bh8s.element_flux_pa_loc_vec[k];  // 0 to 3
//        int b = p_bh8s.element_flux_pb_loc_vec[k];  // 0 to 3
//
//        // identify boundary type
//        int config_id = p_bh8s.element_flux_config_id_vec[k];
//        BoundaryConfigHex8Struct bch8s = p_bh8s.boundary_config_vec[config_id];
//
//        // apply boundary condition
//        if (bch8s.boundary_type_str == "neumann")
//        {
//
//            // calculate distance from point a to b
//            double xa = x_hex8_arr[a];
//            double xb = x_hex8_arr[b];
//            double ya = y_hex8_arr[a];
//            double yb = y_hex8_arr[b];
//            double dist_pa_pb = sqrt((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb));
//
//            // add to b_vec
//            b_vec.coeffRef(start_id_p + n_hex8_arr[a]) += -0.5 * bch8s.boundary_parameter_vec[0] * dist_pa_pb;
//            b_vec.coeffRef(start_id_p + n_hex8_arr[b]) += -0.5 * bch8s.boundary_parameter_vec[0] * dist_pa_pb;
//
//        }
//
//    }

    // clear rows with value boundary elements (x-component of NSE)
    for (int k = 0; k < u_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = u_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = u_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = u_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = u_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = u_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = u_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = u_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = u_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = u_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[c]) = 0.;
        }
        if (d != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[d]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[d]) = 0.;
        }
        if (e != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[e]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[e]) = 0.;
        }
        if (f != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[f]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[f]) = 0.;
        }
        if (g != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[g]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[g]) = 0.;
        }
        if (h != -1)
        {
            a_mat.row(start_id_u + n_hex20_arr[h]) *= 0.;
            b_vec.coeffRef(start_id_u + n_hex20_arr[h]) = 0.;
        }

    }

    // clear rows with value boundary elements (y-component of NSE)
    for (int k = 0; k < v_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = v_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = v_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = v_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = v_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = v_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = v_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = v_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = v_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = v_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[c]) = 0.;
        }
        if (d != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[d]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[d]) = 0.;
        }
        if (e != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[e]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[e]) = 0.;
        }
        if (f != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[f]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[f]) = 0.;
        }
        if (g != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[g]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[g]) = 0.;
        }
        if (h != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[h]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[h]) = 0.;
        }

    }

    // clear rows with value boundary elements (z-component of NSE)
    for (int k = 0; k < w_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = w_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = w_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = w_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = w_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = w_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = w_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = w_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = w_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = w_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_v + n_hex20_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_v + n_hex20_arr[c]) = 0.;
        }
        if (d != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[d]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[d]) = 0.;
        }
        if (e != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[e]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[e]) = 0.;
        }
        if (f != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[f]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[f]) = 0.;
        }
        if (g != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[g]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[g]) = 0.;
        }
        if (h != -1)
        {
            a_mat.row(start_id_w + n_hex20_arr[h]) *= 0.;
            b_vec.coeffRef(start_id_w + n_hex20_arr[h]) = 0.;
        }

    }

    // clear rows with value boundary elements (continuity)
    for (int k = 0; k < p_bh8s.num_element_value; k++)
    {

        // get id of element
        int m = p_bh8s.element_value_id_vec[k];

        // get id of points around hex8 element
        int n00_hex8 = gh8s.element_p00_id_vec[m]; int n04_hex8 = gh8s.element_p04_id_vec[m];
        int n01_hex8 = gh8s.element_p01_id_vec[m]; int n05_hex8 = gh8s.element_p05_id_vec[m];
        int n02_hex8 = gh8s.element_p02_id_vec[m]; int n06_hex8 = gh8s.element_p06_id_vec[m];
        int n03_hex8 = gh8s.element_p03_id_vec[m]; int n07_hex8 = gh8s.element_p07_id_vec[m];
        int n_hex8_arr[8] = {n00_hex8, n01_hex8, n02_hex8, n03_hex8, n04_hex8, n05_hex8, n06_hex8, n07_hex8};

        // get points where the boundary is applied
        int a = p_bh8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = p_bh8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = p_bh8s.element_value_pc_loc_vec[k];  // 0 to 7
        int d = p_bh8s.element_value_pd_loc_vec[k];  // 0 to 7

        // erase entire row
        // -1 values indicate invalid points
        if (a != -1)
        {
            a_mat.row(start_id_p + n_hex8_arr[a]) *= 0.;
            b_vec.coeffRef(start_id_p + n_hex8_arr[a]) = 0.;
        }
        if (b != -1)
        {
            a_mat.row(start_id_p + n_hex8_arr[b]) *= 0.;
            b_vec.coeffRef(start_id_p + n_hex8_arr[b]) = 0.;
        }
        if (c != -1)
        {
            a_mat.row(start_id_p + n_hex8_arr[c]) *= 0.;
            b_vec.coeffRef(start_id_p + n_hex8_arr[c]) = 0.;
        }
        if (d != -1)
        {
            a_mat.row(start_id_p + n_hex8_arr[d]) *= 0.;
            b_vec.coeffRef(start_id_p + n_hex8_arr[d]) = 0.;
        }

    }

    // iterate for each value boundary element (x-component of NSE)
    for (int k = 0; k < u_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = u_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = u_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = u_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = u_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = u_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = u_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = u_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = u_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = u_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // identify boundary type
        int config_id = u_bh20s.element_value_config_id_vec[k];
        BoundaryConfigHex20Struct bch20s = u_bh20s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch20s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[a], start_id_u + n_hex20_arr[a]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[a]) += bch20s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[b], start_id_u + n_hex20_arr[b]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[b]) += bch20s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[c], start_id_u + n_hex20_arr[c]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[c]) += bch20s.boundary_parameter_vec[0];
            }
            if (d != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[d], start_id_u + n_hex20_arr[d]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[d]) += bch20s.boundary_parameter_vec[0];
            }
            if (e != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[e], start_id_u + n_hex20_arr[e]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[e]) += bch20s.boundary_parameter_vec[0];
            }
            if (f != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[f], start_id_u + n_hex20_arr[f]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[f]) += bch20s.boundary_parameter_vec[0];
            }
            if (g != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[g], start_id_u + n_hex20_arr[g]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[g]) += bch20s.boundary_parameter_vec[0];
            }
            if (h != -1)
            {
                a_mat.coeffRef(start_id_u + n_hex20_arr[h], start_id_u + n_hex20_arr[h]) += 1.;
                b_vec.coeffRef(start_id_u + n_hex20_arr[h]) += bch20s.boundary_parameter_vec[0];
            }

        }

    }

    // iterate for each value boundary element (y-component of NSE)
    for (int k = 0; k < v_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = v_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = v_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = v_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = v_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = v_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = v_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = v_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = v_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = v_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // identify boundary type
        int config_id = v_bh20s.element_value_config_id_vec[k];
        BoundaryConfigHex20Struct bch20s = v_bh20s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch20s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[a], start_id_v + n_hex20_arr[a]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[a]) += bch20s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[b], start_id_v + n_hex20_arr[b]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[b]) += bch20s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[c], start_id_v + n_hex20_arr[c]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[c]) += bch20s.boundary_parameter_vec[0];
            }
            if (d != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[d], start_id_v + n_hex20_arr[d]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[d]) += bch20s.boundary_parameter_vec[0];
            }
            if (e != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[e], start_id_v + n_hex20_arr[e]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[e]) += bch20s.boundary_parameter_vec[0];
            }
            if (f != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[f], start_id_v + n_hex20_arr[f]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[f]) += bch20s.boundary_parameter_vec[0];
            }
            if (g != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[g], start_id_v + n_hex20_arr[g]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[g]) += bch20s.boundary_parameter_vec[0];
            }
            if (h != -1)
            {
                a_mat.coeffRef(start_id_v + n_hex20_arr[h], start_id_v + n_hex20_arr[h]) += 1.;
                b_vec.coeffRef(start_id_v + n_hex20_arr[h]) += bch20s.boundary_parameter_vec[0];
            }

        }

    }

    // iterate for each value boundary element (z-component of NSE)
    for (int k = 0; k < w_bh20s.num_element_value; k++)
    {

        // get id of element
        int m = w_bh20s.element_value_id_vec[k];

        // get id of points around hex20 element
        int n00_hex20 = gh20s.element_p00_id_vec[m]; int n10_hex20 = gh20s.element_p10_id_vec[m];
        int n01_hex20 = gh20s.element_p01_id_vec[m]; int n11_hex20 = gh20s.element_p11_id_vec[m];
        int n02_hex20 = gh20s.element_p02_id_vec[m]; int n12_hex20 = gh20s.element_p12_id_vec[m];
        int n03_hex20 = gh20s.element_p03_id_vec[m]; int n13_hex20 = gh20s.element_p13_id_vec[m];
        int n04_hex20 = gh20s.element_p04_id_vec[m]; int n14_hex20 = gh20s.element_p14_id_vec[m];
        int n05_hex20 = gh20s.element_p05_id_vec[m]; int n15_hex20 = gh20s.element_p15_id_vec[m];
        int n06_hex20 = gh20s.element_p06_id_vec[m]; int n16_hex20 = gh20s.element_p16_id_vec[m];
        int n07_hex20 = gh20s.element_p07_id_vec[m]; int n17_hex20 = gh20s.element_p17_id_vec[m];
        int n08_hex20 = gh20s.element_p08_id_vec[m]; int n18_hex20 = gh20s.element_p18_id_vec[m];
        int n09_hex20 = gh20s.element_p09_id_vec[m]; int n19_hex20 = gh20s.element_p19_id_vec[m];
        int n_hex20_arr[20] = {
            n00_hex20, n01_hex20, n02_hex20, n03_hex20, n04_hex20, n05_hex20, n06_hex20, n07_hex20, n08_hex20, n09_hex20,
            n10_hex20, n11_hex20, n12_hex20, n13_hex20, n14_hex20, n15_hex20, n16_hex20, n17_hex20, n18_hex20, n19_hex20
        };

        // get points where the boundary is applied
        int a = w_bh20s.element_value_pa_loc_vec[k];  // 0 to 19
        int b = w_bh20s.element_value_pb_loc_vec[k];  // 0 to 19
        int c = w_bh20s.element_value_pc_loc_vec[k];  // 0 to 19
        int d = w_bh20s.element_value_pd_loc_vec[k];  // 0 to 19
        int e = w_bh20s.element_value_pe_loc_vec[k];  // 0 to 19
        int f = w_bh20s.element_value_pf_loc_vec[k];  // 0 to 19
        int g = w_bh20s.element_value_pg_loc_vec[k];  // 0 to 19
        int h = w_bh20s.element_value_ph_loc_vec[k];  // 0 to 19

        // identify boundary type
        int config_id = w_bh20s.element_value_config_id_vec[k];
        BoundaryConfigHex20Struct bch20s = w_bh20s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch20s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[a], start_id_w + n_hex20_arr[a]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[a]) += bch20s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[b], start_id_w + n_hex20_arr[b]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[b]) += bch20s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[c], start_id_w + n_hex20_arr[c]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[c]) += bch20s.boundary_parameter_vec[0];
            }
            if (d != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[d], start_id_w + n_hex20_arr[d]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[d]) += bch20s.boundary_parameter_vec[0];
            }
            if (e != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[e], start_id_w + n_hex20_arr[e]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[e]) += bch20s.boundary_parameter_vec[0];
            }
            if (f != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[f], start_id_w + n_hex20_arr[f]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[f]) += bch20s.boundary_parameter_vec[0];
            }
            if (g != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[g], start_id_w + n_hex20_arr[g]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[g]) += bch20s.boundary_parameter_vec[0];
            }
            if (h != -1)
            {
                a_mat.coeffRef(start_id_w + n_hex20_arr[h], start_id_w + n_hex20_arr[h]) += 1.;
                b_vec.coeffRef(start_id_w + n_hex20_arr[h]) += bch20s.boundary_parameter_vec[0];
            }

        }

    }

    // iterate for each value boundary element (continuity)
    for (int k = 0; k < p_bh8s.num_element_value; k++)
    {

        // get id of element
        int m = p_bh8s.element_value_id_vec[k];

        // get id of points around hex8 element
        int n00_hex8 = gh8s.element_p00_id_vec[m]; int n04_hex8 = gh8s.element_p04_id_vec[m];
        int n01_hex8 = gh8s.element_p01_id_vec[m]; int n05_hex8 = gh8s.element_p05_id_vec[m];
        int n02_hex8 = gh8s.element_p02_id_vec[m]; int n06_hex8 = gh8s.element_p06_id_vec[m];
        int n03_hex8 = gh8s.element_p03_id_vec[m]; int n07_hex8 = gh8s.element_p07_id_vec[m];
        int n_hex8_arr[8] = {n00_hex8, n01_hex8, n02_hex8, n03_hex8, n04_hex8, n05_hex8, n06_hex8, n07_hex8};

        // get points where the boundary is applied
        int a = p_bh8s.element_value_pa_loc_vec[k];  // 0 to 7
        int b = p_bh8s.element_value_pb_loc_vec[k];  // 0 to 7
        int c = p_bh8s.element_value_pc_loc_vec[k];  // 0 to 7
        int d = p_bh8s.element_value_pd_loc_vec[k];  // 0 to 7

        // identify boundary type
        int config_id = p_bh8s.element_value_config_id_vec[k];
        BoundaryConfigHex8Struct bch8s = p_bh8s.boundary_config_vec[config_id];

        // apply boundary condition
        if (bch8s.boundary_type_str == "dirichlet")
        {

            // set a_mat and b_vec
            // -1 values indicate invalid points
            if (a != -1)
            {
                a_mat.coeffRef(start_id_p + n_hex8_arr[a], start_id_p + n_hex8_arr[a]) += 1.;
                b_vec.coeffRef(start_id_p + n_hex8_arr[a]) += bch8s.boundary_parameter_vec[0];
            }
            if (b != -1)
            {
                a_mat.coeffRef(start_id_p + n_hex8_arr[b], start_id_p + n_hex8_arr[b]) += 1.;
                b_vec.coeffRef(start_id_p + n_hex8_arr[b]) += bch8s.boundary_parameter_vec[0];
            }
            if (c != -1)
            {
                a_mat.coeffRef(start_id_p + n_hex8_arr[c], start_id_p + n_hex8_arr[c]) += 1.;
                b_vec.coeffRef(start_id_p + n_hex8_arr[c]) += bch8s.boundary_parameter_vec[0];
            }
            if (d != -1)
            {
                a_mat.coeffRef(start_id_p + n_hex8_arr[d], start_id_p + n_hex8_arr[d]) += 1.;
                b_vec.coeffRef(start_id_p + n_hex8_arr[d]) += bch8s.boundary_parameter_vec[0];
            }

        }

    }

}

#endif
