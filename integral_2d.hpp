#ifndef INTEGRAL_2D
#define INTEGRAL_2D
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_2d.hpp"

namespace FEM3D
{

class Integral2D
{
    /*

    Test function (N) integrals over a 2D domain.

    Variables
    =========
    domain_in : Domain2D
        Domain where element integrals are calculated.

    Functions
    =========
    evaluate_integral_Ni : void
        Calculates the integral of Ni.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // domain where integral is applied
    Domain2D *domain_ptr;

    // vectors with domain test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_vec;
    VectorDouble3D Ni_vec;
    VectorDouble2D jacobian_determinant_vec;

    // vectors with domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble2D integral_Ni_vec;

    // functions for computing domain integrals
    void evaluate_integral_Ni();

    // default constructor
    Integral2D() {}

    // constructor
    Integral2D(Domain2D &domain_in)
    {
        
        // store domain and boundaries
        domain_ptr = &domain_in;

        // evaluate test functions
        switch (domain_ptr->type_element)
        {
            case 0:
                evaluate_Ni_tri3();
            break;
            case 1:
                evaluate_Ni_quad4();
            break;
        }

    }

    private:
    void evaluate_Ni_tri3();
    void evaluate_Ni_quad4();

};

void Integral2D::evaluate_integral_Ni()
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

void Integral2D::evaluate_Ni_tri3()
{

    // weights for integration
    weight_vec = {0.5};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    const double M_1_3 = 1./3.;
    double a_arr[1] = {M_1_3};
    double b_arr[1] = {M_1_3};

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];

        // get 3d coordinates of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double z0 = domain_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_ptr->point_position_z_vec[p2_pdid];

        // transform 3d points into local 2d coordinates
        Eigen::Vector3d r1_vec;  // p1 relative to p0 in 3d
        Eigen::Vector3d r2_vec;  // p2 relative to p0 in 3d
        r1_vec << x1 - x0, y1 - y0, z1 - z0;
        r2_vec << x2 - x0, y2 - y0, z2 - z0;
        Eigen::Vector3d X_vec;  // unit x coordinate (2d) in 3d coordinates
        Eigen::Vector3d Y_vec;  // unit y coordinate (2d) in 3d coordinates
        Eigen::Vector3d norm_vec;  // outward normal to plane
        X_vec = r2_vec.normalized();
        norm_vec = (r2_vec.cross(r1_vec)).normalized();
        Y_vec = norm_vec.cross(X_vec);
        double X0 = 0;
        double X1 = r1_vec.dot(X_vec);
        double X2 = r2_vec.dot(X_vec);
        double Y0 = 0;
        double Y1 = r1_vec.dot(Y_vec);
        double Y2 = r2_vec.dot(Y_vec);

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 3; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_arr[indx_l];
            double b = b_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_X_a = X0 - X1;
            double derivative_X_b = X2 - X1;
            double derivative_Y_a = Y0 - Y1;
            double derivative_Y_b = Y2 - Y1;

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_X_a, derivative_X_b, derivative_Y_a, derivative_Y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 3; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = a; break;
                    case 1: N = 1. - a - b; break;
                    case 2: N = b; break;
                }

                // store in vectors
                N_part_mli_vec.push_back(N);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_vec.push_back(N_part_ml_vec);
 
    }

}

void Integral2D::evaluate_Ni_quad4()
{

    // weights for integration
    weight_vec = {1., 1., 1., 1.};

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1] * [-1, 1]
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D derivative_N_x_part_ml_vec;
        VectorDouble2D derivative_N_y_part_ml_vec;

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get 3d coordinates of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];
        double z0 = domain_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_ptr->point_position_z_vec[p3_pdid];

        // transform 3d points into local 2d coordinates
        Eigen::Vector3d r1_vec;  // p1 relative to p0 in 3d
        Eigen::Vector3d r2_vec;  // p2 relative to p0 in 3d
        Eigen::Vector3d r3_vec;  // p3 relative to p0 in 3d
        r1_vec << x1 - x0, y1 - y0, z1 - z0;
        r2_vec << x2 - x0, y2 - y0, z2 - z0;
        r3_vec << x3 - x0, y3 - y0, z3 - z0;
        Eigen::Vector3d X_vec;  // unit x coordinate (2d) in 3d coordinates
        Eigen::Vector3d Y_vec;  // unit y coordinate (2d) in 3d coordinates
        Eigen::Vector3d norm_vec;  // outward normal to plane
        X_vec = r3_vec.normalized();
        norm_vec = (r3_vec.cross(r1_vec)).normalized();
        Y_vec = norm_vec.cross(X_vec);
        double X0 = 0;
        double X1 = r1_vec.dot(X_vec);
        double X2 = r2_vec.dot(X_vec);
        double X3 = r3_vec.dot(X_vec);
        double Y0 = 0;
        double Y1 = r1_vec.dot(Y_vec);
        double Y2 = r2_vec.dot(Y_vec);
        double Y3 = r3_vec.dot(Y_vec);

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 4; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_arr[indx_l];
            double b = b_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_X_a = 0.25*(b*X0 - b*X1 + b*X2 - b*X3 - X0 - X1 + X2 + X3);
            double derivative_X_b = 0.25*(a*X0 - a*X1 + a*X2 - a*X3 - X0 + X1 + X2 - X3);
            double derivative_Y_a = 0.25*(b*Y0 - b*Y1 + b*Y2 - b*Y3 - Y0 - Y1 + Y2 + Y3);
            double derivative_Y_b = 0.25*(a*Y0 - a*Y1 + a*Y2 - a*Y3 - Y0 + Y1 + Y2 - Y3);

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_X_a, derivative_X_b, derivative_Y_a, derivative_Y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.25*(1. - a)*(1. - b); break;
                    case 1: N = 0.25*(1. - a)*(1. + b); break;
                    case 2: N = 0.25*(1. + a)*(1. + b); break;
                    case 3: N = 0.25*(1. + a)*(1. - b); break;
                }

                // store in vectors
                N_part_mli_vec.push_back(N);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_vec.push_back(N_part_ml_vec);
 
    }

}

}

#endif
