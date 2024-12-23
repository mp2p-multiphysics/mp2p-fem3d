#ifndef INTEGRAL_3D
#define INTEGRAL_3D
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_3d.hpp"

namespace FEM3D
{

class Integral3D
{
    /*

    Test function (N) integrals over a domain.

    Variables
    =========
    domain_in : Domain3D
        Domain where element integrals are calculated.

    Functions
    =========
    evaluate_integral_Ni : void
        Calculates the integral of Ni.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj.
    evaluate_integral_Ni_derivative_Nj_x : void
        Calculates the integral of Ni * d(Nj)/dx.
    evaluate_integral_Ni_derivative_Nj_y : void
        Calculates the integral of Ni * d(Nj)/dy.
    evaluate_integral_Ni_derivative_Nj_z : void
        Calculates the integral of Ni * d(Nj)/dz.
    evaluate_integral_div_Ni_dot_div_Nj : void
        Calculates the integral of div(Ni) dot div(Nj).

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // domain where integral is applied
    Domain3D *domain_ptr;

    // vectors with domain test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_vec;
    VectorDouble3D Ni_vec;
    VectorDouble3D derivative_Ni_x_vec;
    VectorDouble3D derivative_Ni_y_vec;
    VectorDouble3D derivative_Ni_z_vec;
    VectorDouble2D jacobian_determinant_vec;

    // vectors with domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble2D integral_Ni_vec;
    VectorDouble3D integral_Ni_Nj_vec;
    VectorDouble3D integral_Ni_derivative_Nj_x_vec;
    VectorDouble3D integral_Ni_derivative_Nj_y_vec;
    VectorDouble3D integral_Ni_derivative_Nj_z_vec;
    VectorDouble3D integral_div_Ni_dot_div_Nj_vec;

    // functions for computing domain integrals
    void evaluate_integral_Ni();
    void evaluate_integral_Ni_Nj();
    void evaluate_integral_Ni_derivative_Nj_x();
    void evaluate_integral_Ni_derivative_Nj_y();
    void evaluate_integral_Ni_derivative_Nj_z();
    void evaluate_integral_div_Ni_dot_div_Nj();

    // default constructor
    Integral3D() {}

    // constructor
    Integral3D(Domain3D &domain_in)
    {
        
        // store domain and boundaries
        domain_ptr = &domain_in;

        // evaluate test functions
        switch (domain_ptr->type_element)
        {
            case 0:
                evaluate_Ni_tet4();
            break;
            case 1:
                evaluate_Ni_hex8();
            break;
        }

    }

    private:
    void evaluate_Ni_tet4();
    void evaluate_Ni_hex8();

};

void Integral3D::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Ni_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * Ni_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_integral_Ni_derivative_Nj_x()
{
    /*

    Calculates the integral of Ni * d(Nj)/dx.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_x_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_integral_Ni_derivative_Nj_y()
{
    /*

    Calculates the integral of Ni * d(Nj)/dy.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_y_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_y_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_integral_Ni_derivative_Nj_z()
{
    /*

    Calculates the integral of Ni * d(Nj)/dy.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * Ni_vec[edid][indx_l][indx_i] * derivative_Ni_z_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Nj_z_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_integral_div_Ni_dot_div_Nj()
{
    /*

    Calculates the integral of div(Ni) dot div(Nj).

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_vec.size(); indx_l++) 
        {
            integral_value += weight_vec[indx_l] * jacobian_determinant_vec[edid][indx_l] * (derivative_Ni_x_vec[edid][indx_l][indx_i] * derivative_Ni_x_vec[edid][indx_l][indx_j] + derivative_Ni_y_vec[edid][indx_l][indx_i] * derivative_Ni_y_vec[edid][indx_l][indx_j] + derivative_Ni_z_vec[edid][indx_l][indx_i] * derivative_Ni_z_vec[edid][indx_l][indx_j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Ni_dot_div_Nj_vec.push_back(integral_part_i_vec);

    }

}

void Integral3D::evaluate_Ni_tet4()
{

    // weights for integration
    weight_vec = {0.0416666667, 0.0416666667, 0.0416666667, 0.0416666667};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_arr[4] = {0.5854101966, 0.1381966011, 0.1381966011, 0.1381966011};
    double b_arr[4] = {0.1381966011, 0.5854101966, 0.1381966011, 0.1381966011};
    double c_arr[4] = {0.1381966011, 0.1381966011, 0.5854101966, 0.1381966011};

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];

        // get z values of points
        double z0 = domain_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_ptr->point_position_z_vec[p3_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 4; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;
            VectorDouble derivative_Ni_z_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_arr[indx_l];
            double b = b_arr[indx_l];
            double c = c_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = x1 - x0;
            double derivative_x_b = x2 - x0;
            double derivative_x_c = x3 - x0;
            double derivative_y_a = y1 - y0;
            double derivative_y_b = y2 - y0;
            double derivative_y_c = y3 - y0;
            double derivative_z_a = z1 - z0;
            double derivative_z_b = z2 - z0;
            double derivative_z_c = z3 - z0;

            // get jacobian and its inverse and determinant
            Eigen::Matrix3d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_x_c, derivative_y_a, derivative_y_b, derivative_y_c, derivative_z_a, derivative_z_b, derivative_z_c;
            Eigen::Matrix3d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 1 - a - b - c; break;
                    case 1: N = a; break;
                    case 2: N = b; break;
                    case 3: N = c; break;
                }

                // get derivatives of test function N
                double derivative_Ni_a = 0.;
                double derivative_Ni_b = 0.;
                double derivative_Ni_c = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Ni_a = -1; derivative_Ni_b = -1; derivative_Ni_c = -1; break;
                    case 1: derivative_Ni_a =  1; derivative_Ni_b =  0; derivative_Ni_c =  0; break;
                    case 2: derivative_Ni_a =  0; derivative_Ni_b =  1; derivative_Ni_c =  0; break;
                    case 3: derivative_Ni_a =  0; derivative_Ni_b =  0; derivative_Ni_c =  1; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Ni_abc_vec;
                derivative_Ni_abc_vec << derivative_Ni_a, derivative_Ni_b, derivative_Ni_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Ni_xyz_vec = derivative_Ni_abc_vec * jacobian_inverse_mat;
                double derivative_Ni_x = derivative_Ni_xyz_vec.coeffRef(0);
                double derivative_Ni_y = derivative_Ni_xyz_vec.coeffRef(1);
                double derivative_Ni_z = derivative_Ni_xyz_vec.coeffRef(2);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_Ni_x_part_mli_vec.push_back(derivative_Ni_x);
                derivative_Ni_y_part_mli_vec.push_back(derivative_Ni_y);
                derivative_Ni_z_part_mli_vec.push_back(derivative_Ni_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_vec.push_back(N_part_ml_vec);
        derivative_Ni_x_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_vec.push_back(derivative_Ni_z_part_ml_vec);
 
    }

}

void Integral3D::evaluate_Ni_hex8()
{

    // weights for integration
    weight_vec = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_arr[8] = {-0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692};
    double b_arr[8] = {-0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692};
    double c_arr[8] = {+0.5773502692, +0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, -0.5773502692, -0.5773502692};

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p4_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][4];
        int p5_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][5];
        int p6_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][6];
        int p7_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][7];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];
        int p4_pdid = domain_ptr->point_pgid_to_pdid_map[p4_pgid];
        int p5_pdid = domain_ptr->point_pgid_to_pdid_map[p5_pgid];
        int p6_pdid = domain_ptr->point_pgid_to_pdid_map[p6_pgid];
        int p7_pdid = domain_ptr->point_pgid_to_pdid_map[p7_pgid];

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];
        double x4 = domain_ptr->point_position_x_vec[p4_pdid];
        double x5 = domain_ptr->point_position_x_vec[p5_pdid];
        double x6 = domain_ptr->point_position_x_vec[p6_pdid];
        double x7 = domain_ptr->point_position_x_vec[p7_pdid];

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];
        double y4 = domain_ptr->point_position_y_vec[p4_pdid];
        double y5 = domain_ptr->point_position_y_vec[p5_pdid];
        double y6 = domain_ptr->point_position_y_vec[p6_pdid];
        double y7 = domain_ptr->point_position_y_vec[p7_pdid];

        // get z values of points
        double z0 = domain_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_ptr->point_position_z_vec[p3_pdid];
        double z4 = domain_ptr->point_position_z_vec[p4_pdid];
        double z5 = domain_ptr->point_position_z_vec[p5_pdid];
        double z6 = domain_ptr->point_position_z_vec[p6_pdid];
        double z7 = domain_ptr->point_position_z_vec[p7_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 8; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;
            VectorDouble derivative_Ni_z_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_arr[indx_l];
            double b = b_arr[indx_l];
            double c = c_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = -x0*(b - 1)*(c - 1)*0.125 + x1*(b - 1)*(c - 1)*0.125 - x2*(b + 1)*(c - 1)*0.125 + x3*(b + 1)*(c - 1)*0.125 + x4*(b - 1)*(c + 1)*0.125 - x5*(b - 1)*(c + 1)*0.125 + x6*(b + 1)*(c + 1)*0.125 - x7*(b + 1)*(c + 1)*0.125;
            double derivative_x_b = -x0*(a - 1)*(c - 1)*0.125 + x1*(a + 1)*(c - 1)*0.125 - x2*(a + 1)*(c - 1)*0.125 + x3*(a - 1)*(c - 1)*0.125 + x4*(a - 1)*(c + 1)*0.125 - x5*(a + 1)*(c + 1)*0.125 + x6*(a + 1)*(c + 1)*0.125 - x7*(a - 1)*(c + 1)*0.125;
            double derivative_x_c = -x0*(a - 1)*(b - 1)*0.125 + x1*(a + 1)*(b - 1)*0.125 - x2*(a + 1)*(b + 1)*0.125 + x3*(a - 1)*(b + 1)*0.125 + x4*(a - 1)*(b - 1)*0.125 - x5*(a + 1)*(b - 1)*0.125 + x6*(a + 1)*(b + 1)*0.125 - x7*(a - 1)*(b + 1)*0.125;
            double derivative_y_a = -y0*(b - 1)*(c - 1)*0.125 + y1*(b - 1)*(c - 1)*0.125 - y2*(b + 1)*(c - 1)*0.125 + y3*(b + 1)*(c - 1)*0.125 + y4*(b - 1)*(c + 1)*0.125 - y5*(b - 1)*(c + 1)*0.125 + y6*(b + 1)*(c + 1)*0.125 - y7*(b + 1)*(c + 1)*0.125;
            double derivative_y_b = -y0*(a - 1)*(c - 1)*0.125 + y1*(a + 1)*(c - 1)*0.125 - y2*(a + 1)*(c - 1)*0.125 + y3*(a - 1)*(c - 1)*0.125 + y4*(a - 1)*(c + 1)*0.125 - y5*(a + 1)*(c + 1)*0.125 + y6*(a + 1)*(c + 1)*0.125 - y7*(a - 1)*(c + 1)*0.125;
            double derivative_y_c = -y0*(a - 1)*(b - 1)*0.125 + y1*(a + 1)*(b - 1)*0.125 - y2*(a + 1)*(b + 1)*0.125 + y3*(a - 1)*(b + 1)*0.125 + y4*(a - 1)*(b - 1)*0.125 - y5*(a + 1)*(b - 1)*0.125 + y6*(a + 1)*(b + 1)*0.125 - y7*(a - 1)*(b + 1)*0.125;
            double derivative_z_a = -z0*(b - 1)*(c - 1)*0.125 + z1*(b - 1)*(c - 1)*0.125 - z2*(b + 1)*(c - 1)*0.125 + z3*(b + 1)*(c - 1)*0.125 + z4*(b - 1)*(c + 1)*0.125 - z5*(b - 1)*(c + 1)*0.125 + z6*(b + 1)*(c + 1)*0.125 - z7*(b + 1)*(c + 1)*0.125;
            double derivative_z_b = -z0*(a - 1)*(c - 1)*0.125 + z1*(a + 1)*(c - 1)*0.125 - z2*(a + 1)*(c - 1)*0.125 + z3*(a - 1)*(c - 1)*0.125 + z4*(a - 1)*(c + 1)*0.125 - z5*(a + 1)*(c + 1)*0.125 + z6*(a + 1)*(c + 1)*0.125 - z7*(a - 1)*(c + 1)*0.125;
            double derivative_z_c = -z0*(a - 1)*(b - 1)*0.125 + z1*(a + 1)*(b - 1)*0.125 - z2*(a + 1)*(b + 1)*0.125 + z3*(a - 1)*(b + 1)*0.125 + z4*(a - 1)*(b - 1)*0.125 - z5*(a + 1)*(b - 1)*0.125 + z6*(a + 1)*(b + 1)*0.125 - z7*(a - 1)*(b + 1)*0.125;

            // get jacobian and its inverse and determinant
            Eigen::Matrix3d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_x_c, derivative_y_a, derivative_y_b, derivative_y_c, derivative_z_a, derivative_z_b, derivative_z_c;
            Eigen::Matrix3d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 8; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = -(a - 1)*(b - 1)*(c - 1)*0.125; break;
                    case 1: N =  (a + 1)*(b - 1)*(c - 1)*0.125; break;
                    case 2: N = -(a + 1)*(b + 1)*(c - 1)*0.125; break;
                    case 3: N =  (a - 1)*(b + 1)*(c - 1)*0.125; break;
                    case 4: N =  (a - 1)*(b - 1)*(c + 1)*0.125; break;
                    case 5: N = -(a + 1)*(b - 1)*(c + 1)*0.125; break;
                    case 6: N =  (a + 1)*(b + 1)*(c + 1)*0.125; break;
                    case 7: N = -(a - 1)*(b + 1)*(c + 1)*0.125; break;
                }

                // get derivatives of test function N
                double derivative_Ni_a = 0.;
                double derivative_Ni_b = 0.;
                double derivative_Ni_c = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Ni_a = -(b - 1)*(c - 1)*0.125; derivative_Ni_b = -(a - 1)*(c - 1)*0.125; derivative_Ni_c = -(a - 1)*(b - 1)*0.125; break;
                    case 1: derivative_Ni_a =  (b - 1)*(c - 1)*0.125; derivative_Ni_b =  (a + 1)*(c - 1)*0.125; derivative_Ni_c =  (a + 1)*(b - 1)*0.125; break;
                    case 2: derivative_Ni_a = -(b + 1)*(c - 1)*0.125; derivative_Ni_b = -(a + 1)*(c - 1)*0.125; derivative_Ni_c = -(a + 1)*(b + 1)*0.125; break;
                    case 3: derivative_Ni_a =  (b + 1)*(c - 1)*0.125; derivative_Ni_b =  (a - 1)*(c - 1)*0.125; derivative_Ni_c =  (a - 1)*(b + 1)*0.125; break;
                    case 4: derivative_Ni_a =  (b - 1)*(c + 1)*0.125; derivative_Ni_b =  (a - 1)*(c + 1)*0.125; derivative_Ni_c =  (a - 1)*(b - 1)*0.125; break;
                    case 5: derivative_Ni_a = -(b - 1)*(c + 1)*0.125; derivative_Ni_b = -(a + 1)*(c + 1)*0.125; derivative_Ni_c = -(a + 1)*(b - 1)*0.125; break;
                    case 6: derivative_Ni_a =  (b + 1)*(c + 1)*0.125; derivative_Ni_b =  (a + 1)*(c + 1)*0.125; derivative_Ni_c =  (a + 1)*(b + 1)*0.125; break;
                    case 7: derivative_Ni_a = -(b + 1)*(c + 1)*0.125; derivative_Ni_b = -(a - 1)*(c + 1)*0.125; derivative_Ni_c = -(a - 1)*(b + 1)*0.125; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Ni_abc_vec;
                derivative_Ni_abc_vec << derivative_Ni_a, derivative_Ni_b, derivative_Ni_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Ni_xyz_vec = derivative_Ni_abc_vec * jacobian_inverse_mat;
                double derivative_Ni_x = derivative_Ni_xyz_vec.coeffRef(0);
                double derivative_Ni_y = derivative_Ni_xyz_vec.coeffRef(1);
                double derivative_Ni_z = derivative_Ni_xyz_vec.coeffRef(2);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_Ni_x_part_mli_vec.push_back(derivative_Ni_x);
                derivative_Ni_y_part_mli_vec.push_back(derivative_Ni_y);
                derivative_Ni_z_part_mli_vec.push_back(derivative_Ni_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_vec.push_back(N_part_ml_vec);
        derivative_Ni_x_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_vec.push_back(derivative_Ni_z_part_ml_vec);
 
    }

}

}

#endif
