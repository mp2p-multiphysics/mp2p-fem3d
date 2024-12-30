#ifndef INTEGRAL_TAYLORHOOD_3D
#define INTEGRAL_TAYLORHOOD_3D
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_3d.hpp"

namespace FEM3D
{

class IntegralTaylorHood3D
{
    /*

    Taylor-Hood test function (mixed linear (N) and quadratic (M)) integrals over a domain.

    Variables
    =========
    domain_line_in : Domain3D
        Linear domain where element integrals are calculated.
    domain_quad_in : Domain3D
        Quadratic domain where element integrals are calculated.

    Functions
    =========
    evaluate_integral_Ni_derivative_Mj_x_line : void
        Calculates the integral of Ni * d(Mj)/dx over test functions in the linear domain.
    evaluate_integral_Ni_derivative_Mj_y_line : void
        Calculates the integral of Ni * d(Mj)/dy over test functions in the linear domain.
    evaluate_integral_Ni_derivative_Mj_z_line : void
        Calculates the integral of Ni * d(Mj)/dz over test functions in the linear domain.
    evaluate_integral_Mi_quad : void
        Calculates the integral of Mi over test functions in the quadratic domain.
    evaluate_integral_Mi_Mj_quad : void
        Calculates the integral of Mi * Mj over test functions in the quadratic domain.
    evaluate_integral_Mi_derivative_Nj_x_quad : void
        Calculates the integral of Mi * d(Nj)/dx over test functions in the quadratic domain.
    evaluate_integral_Mi_derivative_Nj_y_quad : void
        Calculates the integral of Mi * d(Nj)/dy over test functions in the quadratic domain.
    evaluate_integral_Mi_derivative_Nj_z_quad : void
        Calculates the integral of Mi * d(Nj)/dz over test functions in the quadratic domain.
    evaluate_integral_div_Mi_dot_div_Mj_quad : void
        Calculates the integral of div(Mi) dot div(Mj) over test functions in the quadratic domain.
    evaluate_integral_Mi_Mk_derivative_Mj_x_quad : void
        Calculates the integral of Mi * Mk * d(Mj)/dx over test functions in the quadratic domain.
    evaluate_integral_Mi_Mk_derivative_Mj_y_quad : void
        Calculates the integral of Mi * Mk * d(Mj)/dy over test functions in the quadratic domain.
    evaluate_integral_Mi_Mk_derivative_Mj_z_quad : void
        Calculates the integral of Mi * Mk * d(Mj)/dz over test functions in the quadratic domain.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // domain where integral is applied
    Domain3D *domain_line_ptr;
    Domain3D *domain_quad_ptr;

    // vectors with linear test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_line_vec;
    VectorDouble3D Ni_line_vec;
    VectorDouble3D Mi_line_vec;
    VectorDouble3D derivative_Ni_x_line_vec;
    VectorDouble3D derivative_Ni_y_line_vec;
    VectorDouble3D derivative_Ni_z_line_vec;
    VectorDouble3D derivative_Mi_x_line_vec;
    VectorDouble3D derivative_Mi_y_line_vec;
    VectorDouble3D derivative_Mi_z_line_vec;
    VectorDouble2D jacobian_determinant_line_vec;

    // vectors with quadratic test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_quad_vec;
    VectorDouble3D Ni_quad_vec;
    VectorDouble3D Mi_quad_vec;
    VectorDouble3D derivative_Ni_x_quad_vec;
    VectorDouble3D derivative_Ni_y_quad_vec;
    VectorDouble3D derivative_Ni_z_quad_vec;
    VectorDouble3D derivative_Mi_x_quad_vec;
    VectorDouble3D derivative_Mi_y_quad_vec;
    VectorDouble3D derivative_Mi_z_quad_vec;
    VectorDouble2D jacobian_determinant_quad_vec;

    // vectors with linear domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble3D integral_Ni_derivative_Mj_x_line_vec;
    VectorDouble3D integral_Ni_derivative_Mj_y_line_vec;
    VectorDouble3D integral_Ni_derivative_Mj_z_line_vec;
    VectorDouble3D integral_Ni_Nj_line_vec;

    // vectors with quadratic domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble2D integral_Mi_quad_vec;
    VectorDouble3D integral_Mi_Mj_quad_vec;
    VectorDouble3D integral_Mi_derivative_Nj_x_quad_vec;
    VectorDouble3D integral_Mi_derivative_Nj_y_quad_vec;
    VectorDouble3D integral_Mi_derivative_Nj_z_quad_vec;
    VectorDouble3D integral_div_Mi_dot_div_Mj_quad_vec;
    VectorDouble4D integral_Mi_Mk_derivative_Mj_x_quad_vec;
    VectorDouble4D integral_Mi_Mk_derivative_Mj_y_quad_vec;
    VectorDouble4D integral_Mi_Mk_derivative_Mj_z_quad_vec;

    // functions for computing linear domain integrals
    void evaluate_integral_Ni_derivative_Mj_x_line();
    void evaluate_integral_Ni_derivative_Mj_y_line();
    void evaluate_integral_Ni_derivative_Mj_z_line();
    void evaluate_integral_Ni_Nj_line();

    // functions for computing quadratic domain integrals
    void evaluate_integral_Mi_quad();
    void evaluate_integral_Mi_Mj_quad();
    void evaluate_integral_Mi_derivative_Nj_x_quad();
    void evaluate_integral_Mi_derivative_Nj_y_quad();
    void evaluate_integral_Mi_derivative_Nj_z_quad();
    void evaluate_integral_div_Mi_dot_div_Mj_quad();
    void evaluate_integral_Mi_Mk_derivative_Mj_x_quad();
    void evaluate_integral_Mi_Mk_derivative_Mj_y_quad();
    void evaluate_integral_Mi_Mk_derivative_Mj_z_quad();

    // default constructor
    IntegralTaylorHood3D() {}

    // constructor
    IntegralTaylorHood3D(Domain3D &domain_line_in, Domain3D &domain_quad_in)
    {
        
        // store domain and boundaries
        domain_line_ptr = &domain_line_in;
        domain_quad_ptr = &domain_quad_in;

        // evaluate test functions
        if (domain_line_ptr->type_element == 0 && domain_quad_ptr->type_element == 2)
        {
            evaluate_Mi_Ni_tet();
        }
        else if (domain_line_ptr->type_element == 1 && domain_quad_ptr->type_element == 3)
        {
            evaluate_Mi_Ni_hex();
        }

    }

    private:
    void evaluate_Mi_Ni_tet();
    void evaluate_Mi_Ni_hex();

};

void IntegralTaylorHood3D::evaluate_integral_Ni_derivative_Mj_x_line()
{
    /*

    Calculates the integral of Ni * d(Mj)/dx over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * derivative_Mi_x_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Mj_x_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Ni_derivative_Mj_y_line()
{
    /*

    Calculates the integral of Ni * d(Mj)/dy over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * derivative_Mi_y_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Mj_y_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Ni_derivative_Mj_z_line()
{
    /*

    Calculates the integral of Ni * d(Mj)/dz over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * derivative_Mi_z_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Mj_z_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Ni_Nj_line()
{
    /*

    Calculates the integral of Ni * Nj over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * Ni_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_quad()
{
    /*

    Calculates the integral of Mi over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Mi_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_Mj_quad()
{
    /*

    Calculates the integral of Mi * Mj over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mj_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_derivative_Nj_x_quad()
{
    /*

    Calculates the integral of Mi * d(Nj)/dx over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * derivative_Ni_x_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_derivative_Nj_x_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_derivative_Nj_y_quad()
{
    /*

    Calculates the integral of Mi * d(Nj)/dy over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * derivative_Ni_y_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_derivative_Nj_y_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_derivative_Nj_z_quad()
{
    /*

    Calculates the integral of Mi * d(Nj)/dz over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * derivative_Ni_z_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_derivative_Nj_z_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_div_Mi_dot_div_Mj_quad()
{
    /*

    Calculates the integral of div(Mi) dot div(Mj) over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {  
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * (derivative_Mi_x_quad_vec[edid][indx_l][indx_i] * derivative_Mi_x_quad_vec[edid][indx_l][indx_j] + derivative_Mi_y_quad_vec[edid][indx_l][indx_i] * derivative_Mi_y_quad_vec[edid][indx_l][indx_j] + derivative_Mi_z_quad_vec[edid][indx_l][indx_i] * derivative_Mi_z_quad_vec[edid][indx_l][indx_j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Mi_dot_div_Mj_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_Mk_derivative_Mj_x_quad()
{
    /*

    Calculates the integral of Mi * Mk * d(Mj)/dx over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */ 

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_k] * derivative_Mi_x_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mk_derivative_Mj_x_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_Mk_derivative_Mj_y_quad()
{
    /*

    Calculates the integral of Mi * Mk * d(Mj)/dy over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */ 

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_k] * derivative_Mi_y_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mk_derivative_Mj_y_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_integral_Mi_Mk_derivative_Mj_z_quad()
{
    /*

    Calculates the integral of Mi * Mk * d(Mj)/dz over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */ 

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_k] * derivative_Mi_z_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mk_derivative_Mj_z_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood3D::evaluate_Mi_Ni_tet()
{

    // linear integration points

    // weights for integration
    weight_line_vec = {0.0416666667, 0.0416666667, 0.0416666667, 0.0416666667};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_arr[4] = {0.5854101966, 0.1381966011, 0.1381966011, 0.1381966011};
    double b_arr[4] = {0.1381966011, 0.5854101966, 0.1381966011, 0.1381966011};
    double c_arr[4] = {0.1381966011, 0.1381966011, 0.5854101966, 0.1381966011};

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;
        VectorDouble2D derivative_Mi_z_part_ml_vec;

        // get points around element
        int p0_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p0_pdid = domain_line_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_line_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_line_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_line_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get x values of points
        double x0 = domain_line_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_line_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_line_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_line_ptr->point_position_x_vec[p3_pdid];

        // get y values of points
        double y0 = domain_line_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_line_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_line_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_line_ptr->point_position_y_vec[p3_pdid];

        // get z values of points
        double z0 = domain_line_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_line_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_line_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_line_ptr->point_position_z_vec[p3_pdid];
 
        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 4; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Mi_z_part_mli_vec;
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

            // iterate for each test function
            for (int indx_i = 0; indx_i < 10; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = (a + b + c - 1)*(2*a + 2*b + 2*c - 1); break;
                    case 1: M = a*(2*a - 1); break;
                    case 2: M = b*(2*b - 1); break;
                    case 3: M = c*(2*c - 1); break;
                    case 4: M = 4*a*(-a - b - c + 1); break;
                    case 5: M = 4*b*c; break;
                    case 6: M = 4*b*(-a - b - c + 1); break;
                    case 7: M = 4*c*(-a - b - c + 1); break;
                    case 8: M = 4*a*b; break;
                    case 9: M = 4*a*c; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                double derivative_Mi_c = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Mi_a = 4*a + 4*b + 4*c - 3;  derivative_Mi_b = 4*a + 4*b + 4*c - 3;  derivative_Mi_c = 4*a + 4*b + 4*c - 3;  break;
                    case 1: derivative_Mi_a = 4*a - 1;              derivative_Mi_b = 0;                    derivative_Mi_c = 0;                    break;
                    case 2: derivative_Mi_a = 0;                    derivative_Mi_b = 4*b - 1;              derivative_Mi_c = 0;                    break;
                    case 3: derivative_Mi_a = 0;                    derivative_Mi_b = 0;                    derivative_Mi_c = 4*c - 1;              break;
                    case 4: derivative_Mi_a = -8*a - 4*b - 4*c + 4; derivative_Mi_b = -4*a;                 derivative_Mi_c = -4*a;                 break;
                    case 5: derivative_Mi_a = 0;                    derivative_Mi_b = 4*c;                  derivative_Mi_c = 4*b;                  break;
                    case 6: derivative_Mi_a = -4*b;                 derivative_Mi_b = -4*a - 8*b - 4*c + 4; derivative_Mi_c = -4*b;                 break;
                    case 7: derivative_Mi_a = -4*c;                 derivative_Mi_b = -4*c;                 derivative_Mi_c = -4*a - 4*b - 8*c + 4; break;
                    case 8: derivative_Mi_a = 4*b;                  derivative_Mi_b = 4*a;                  derivative_Mi_c = 0;                    break;
                    case 9: derivative_Mi_a = 4*c;                  derivative_Mi_b = 0;                    derivative_Mi_c = 4*a;                  break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Mi_abc_vec;
                derivative_Mi_abc_vec << derivative_Mi_a, derivative_Mi_b, derivative_Mi_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Mi_xyz_vec = derivative_Mi_abc_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xyz_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xyz_vec.coeffRef(1);
                double derivative_Mi_z = derivative_Mi_xyz_vec.coeffRef(2);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);
                derivative_Mi_z_part_mli_vec.push_back(derivative_Mi_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);
            derivative_Mi_z_part_ml_vec.push_back(derivative_Mi_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_line_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_line_vec.push_back(N_part_ml_vec);
        Mi_line_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_line_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_line_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_line_vec.push_back(derivative_Ni_z_part_ml_vec);
        derivative_Mi_x_line_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_line_vec.push_back(derivative_Mi_y_part_ml_vec);
        derivative_Mi_z_line_vec.push_back(derivative_Mi_z_part_ml_vec);
 
    }

    // quadratic integration points

    // weights for integration
    weight_quad_vec = {0.0027777778, 0.0027777778, 0.0027777778, 0.0027777778, 0.0111111111, 0.0111111111, 0.0111111111, 0.0111111111, 0.0111111111, 0.0111111111, 0.0888888889};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_quad_arr[11] = {1.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.50, 0.00, 0.50, 0.25};
    double b_quad_arr[11] = {0.00, 1.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.50, 0.00, 0.25};
    double c_quad_arr[11] = {0.00, 0.00, 1.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.50, 0.25};

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;
        VectorDouble2D derivative_Mi_z_part_ml_vec;

        // get points around element
        int p0_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p4_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][4];
        int p5_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][5];
        int p6_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][6];
        int p7_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][7];
        int p8_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][8];
        int p9_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][9];
        int p0_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p3_pgid];
        int p4_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p4_pgid];
        int p5_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p5_pgid];
        int p6_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p6_pgid];
        int p7_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p7_pgid];
        int p8_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p8_pgid];
        int p9_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p9_pgid];

        // get x values of points
        double x0 = domain_quad_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_quad_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_quad_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_quad_ptr->point_position_x_vec[p3_pdid];
        double x4 = domain_quad_ptr->point_position_x_vec[p4_pdid];
        double x5 = domain_quad_ptr->point_position_x_vec[p5_pdid];
        double x6 = domain_quad_ptr->point_position_x_vec[p6_pdid];
        double x7 = domain_quad_ptr->point_position_x_vec[p7_pdid];
        double x8 = domain_quad_ptr->point_position_x_vec[p8_pdid];
        double x9 = domain_quad_ptr->point_position_x_vec[p9_pdid];

        // get y values of points
        double y0 = domain_quad_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_quad_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_quad_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_quad_ptr->point_position_y_vec[p3_pdid];
        double y4 = domain_quad_ptr->point_position_y_vec[p4_pdid];
        double y5 = domain_quad_ptr->point_position_y_vec[p5_pdid];
        double y6 = domain_quad_ptr->point_position_y_vec[p6_pdid];
        double y7 = domain_quad_ptr->point_position_y_vec[p7_pdid];
        double y8 = domain_quad_ptr->point_position_y_vec[p8_pdid];
        double y9 = domain_quad_ptr->point_position_y_vec[p9_pdid];

        // get z values of points
        double z0 = domain_quad_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_quad_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_quad_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_quad_ptr->point_position_z_vec[p3_pdid];
        double z4 = domain_quad_ptr->point_position_z_vec[p4_pdid];
        double z5 = domain_quad_ptr->point_position_z_vec[p5_pdid];
        double z6 = domain_quad_ptr->point_position_z_vec[p6_pdid];
        double z7 = domain_quad_ptr->point_position_z_vec[p7_pdid];
        double z8 = domain_quad_ptr->point_position_z_vec[p8_pdid];
        double z9 = domain_quad_ptr->point_position_z_vec[p9_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 11; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Mi_z_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;
            VectorDouble derivative_Ni_z_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_quad_arr[indx_l];
            double b = b_quad_arr[indx_l];
            double c = c_quad_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = 4*a*x0 + 4*a*x1 - 8*a*x4 + 4*b*x0 - 4*b*x4 - 4*b*x6 + 4*b*x8 + 4*c*x0 - 4*c*x4 - 4*c*x7 + 4*c*x9 - 3*x0 - x1 + 4*x4;
            double derivative_x_b = 4*a*x0 - 4*a*x4 - 4*a*x6 + 4*a*x8 + 4*b*x0 + 4*b*x2 - 8*b*x6 + 4*c*x0 + 4*c*x5 - 4*c*x6 - 4*c*x7 - 3*x0 - x2 + 4*x6;
            double derivative_x_c = 4*a*x0 - 4*a*x4 - 4*a*x7 + 4*a*x9 + 4*b*x0 + 4*b*x5 - 4*b*x6 - 4*b*x7 + 4*c*x0 + 4*c*x3 - 8*c*x7 - 3*x0 - x3 + 4*x7;
            double derivative_y_a = 4*a*y0 + 4*a*y1 - 8*a*y4 + 4*b*y0 - 4*b*y4 - 4*b*y6 + 4*b*y8 + 4*c*y0 - 4*c*y4 - 4*c*y7 + 4*c*y9 - 3*y0 - y1 + 4*y4;
            double derivative_y_b = 4*a*y0 - 4*a*y4 - 4*a*y6 + 4*a*y8 + 4*b*y0 + 4*b*y2 - 8*b*y6 + 4*c*y0 + 4*c*y5 - 4*c*y6 - 4*c*y7 - 3*y0 - y2 + 4*y6;
            double derivative_y_c = 4*a*y0 - 4*a*y4 - 4*a*y7 + 4*a*y9 + 4*b*y0 + 4*b*y5 - 4*b*y6 - 4*b*y7 + 4*c*y0 + 4*c*y3 - 8*c*y7 - 3*y0 - y3 + 4*y7;
            double derivative_z_a = 4*a*z0 + 4*a*z1 - 8*a*z4 + 4*b*z0 - 4*b*z4 - 4*b*z6 + 4*b*z8 + 4*c*z0 - 4*c*z4 - 4*c*z7 + 4*c*z9 - 3*z0 - z1 + 4*z4;
            double derivative_z_b = 4*a*z0 - 4*a*z4 - 4*a*z6 + 4*a*z8 + 4*b*z0 + 4*b*z2 - 8*b*z6 + 4*c*z0 + 4*c*z5 - 4*c*z6 - 4*c*z7 - 3*z0 - z2 + 4*z6;
            double derivative_z_c = 4*a*z0 - 4*a*z4 - 4*a*z7 + 4*a*z9 + 4*b*z0 + 4*b*z5 - 4*b*z6 - 4*b*z7 + 4*c*z0 + 4*c*z3 - 8*c*z7 - 3*z0 - z3 + 4*z7;

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

            // iterate for each test function
            for (int indx_i = 0; indx_i < 10; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = (a + b + c - 1)*(2*a + 2*b + 2*c - 1); break;
                    case 1: M = a*(2*a - 1); break;
                    case 2: M = b*(2*b - 1); break;
                    case 3: M = c*(2*c - 1); break;
                    case 4: M = 4*a*(-a - b - c + 1); break;
                    case 5: M = 4*b*c; break;
                    case 6: M = 4*b*(-a - b - c + 1); break;
                    case 7: M = 4*c*(-a - b - c + 1); break;
                    case 8: M = 4*a*b; break;
                    case 9: M = 4*a*c; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                double derivative_Mi_c = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Mi_a = 4*a + 4*b + 4*c - 3;  derivative_Mi_b = 4*a + 4*b + 4*c - 3;  derivative_Mi_c = 4*a + 4*b + 4*c - 3;  break;
                    case 1: derivative_Mi_a = 4*a - 1;              derivative_Mi_b = 0;                    derivative_Mi_c = 0;                    break;
                    case 2: derivative_Mi_a = 0;                    derivative_Mi_b = 4*b - 1;              derivative_Mi_c = 0;                    break;
                    case 3: derivative_Mi_a = 0;                    derivative_Mi_b = 0;                    derivative_Mi_c = 4*c - 1;              break;
                    case 4: derivative_Mi_a = -8*a - 4*b - 4*c + 4; derivative_Mi_b = -4*a;                 derivative_Mi_c = -4*a;                 break;
                    case 5: derivative_Mi_a = 0;                    derivative_Mi_b = 4*c;                  derivative_Mi_c = 4*b;                  break;
                    case 6: derivative_Mi_a = -4*b;                 derivative_Mi_b = -4*a - 8*b - 4*c + 4; derivative_Mi_c = -4*b;                 break;
                    case 7: derivative_Mi_a = -4*c;                 derivative_Mi_b = -4*c;                 derivative_Mi_c = -4*a - 4*b - 8*c + 4; break;
                    case 8: derivative_Mi_a = 4*b;                  derivative_Mi_b = 4*a;                  derivative_Mi_c = 0;                    break;
                    case 9: derivative_Mi_a = 4*c;                  derivative_Mi_b = 0;                    derivative_Mi_c = 4*a;                  break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Mi_abc_vec;
                derivative_Mi_abc_vec << derivative_Mi_a, derivative_Mi_b, derivative_Mi_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Mi_xyz_vec = derivative_Mi_abc_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xyz_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xyz_vec.coeffRef(1);
                double derivative_Mi_z = derivative_Mi_xyz_vec.coeffRef(2);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);
                derivative_Mi_z_part_mli_vec.push_back(derivative_Mi_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);
            derivative_Mi_z_part_ml_vec.push_back(derivative_Mi_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_quad_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_quad_vec.push_back(N_part_ml_vec);
        Mi_quad_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_quad_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_quad_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_quad_vec.push_back(derivative_Ni_z_part_ml_vec);
        derivative_Mi_x_quad_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_quad_vec.push_back(derivative_Mi_y_part_ml_vec);
        derivative_Mi_z_quad_vec.push_back(derivative_Mi_z_part_ml_vec);
 
    }

}

void IntegralTaylorHood3D::evaluate_Mi_Ni_hex()
{

    // linear integration points

    // weights for integration
    weight_line_vec = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1] * [-1, 1]
    double a_arr[8] = {-0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692};
    double b_arr[8] = {-0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, +0.5773502692, +0.5773502692};
    double c_arr[8] = {+0.5773502692, +0.5773502692, +0.5773502692, +0.5773502692, -0.5773502692, -0.5773502692, -0.5773502692, -0.5773502692};

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;
        VectorDouble2D derivative_Mi_z_part_ml_vec;
        
        // get points around element
        int p0_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p4_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][4];
        int p5_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][5];
        int p6_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][6];
        int p7_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][7];
        int p0_pdid = domain_line_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_line_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_line_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_line_ptr->point_pgid_to_pdid_map[p3_pgid];
        int p4_pdid = domain_line_ptr->point_pgid_to_pdid_map[p4_pgid];
        int p5_pdid = domain_line_ptr->point_pgid_to_pdid_map[p5_pgid];
        int p6_pdid = domain_line_ptr->point_pgid_to_pdid_map[p6_pgid];
        int p7_pdid = domain_line_ptr->point_pgid_to_pdid_map[p7_pgid];

        // get x values of points
        double x0 = domain_line_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_line_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_line_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_line_ptr->point_position_x_vec[p3_pdid];
        double x4 = domain_line_ptr->point_position_x_vec[p4_pdid];
        double x5 = domain_line_ptr->point_position_x_vec[p5_pdid];
        double x6 = domain_line_ptr->point_position_x_vec[p6_pdid];
        double x7 = domain_line_ptr->point_position_x_vec[p7_pdid];

        // get y values of points
        double y0 = domain_line_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_line_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_line_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_line_ptr->point_position_y_vec[p3_pdid];
        double y4 = domain_line_ptr->point_position_y_vec[p4_pdid];
        double y5 = domain_line_ptr->point_position_y_vec[p5_pdid];
        double y6 = domain_line_ptr->point_position_y_vec[p6_pdid];
        double y7 = domain_line_ptr->point_position_y_vec[p7_pdid];

        // get z values of points
        double z0 = domain_line_ptr->point_position_z_vec[p0_pdid];
        double z1 = domain_line_ptr->point_position_z_vec[p1_pdid];
        double z2 = domain_line_ptr->point_position_z_vec[p2_pdid];
        double z3 = domain_line_ptr->point_position_z_vec[p3_pdid];
        double z4 = domain_line_ptr->point_position_z_vec[p4_pdid];
        double z5 = domain_line_ptr->point_position_z_vec[p5_pdid];
        double z6 = domain_line_ptr->point_position_z_vec[p6_pdid];
        double z7 = domain_line_ptr->point_position_z_vec[p7_pdid];
 
        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 8; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Mi_z_part_mli_vec;
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

            // iterate for each test function
            for (int indx_i = 0; indx_i < 20; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case  0: M =  (a - 1)*(b - 1)*(c - 1)*(a + b + c + 2)*0.125; break;
                    case  1: M = -(a + 1)*(b - 1)*(c - 1)*(-a + b + c + 2)*0.125; break;
                    case  2: M = -(a + 1)*(b + 1)*(c - 1)*(a + b - c - 2)*0.125; break;
                    case  3: M = -(a - 1)*(b + 1)*(c - 1)*(a - b + c + 2)*0.125; break;
                    case  4: M = -(a - 1)*(b - 1)*(c + 1)*(a + b - c + 2)*0.125; break;
                    case  5: M = -(a + 1)*(b - 1)*(c + 1)*(a - b + c - 2)*0.125; break;
                    case  6: M =  (a + 1)*(b + 1)*(c + 1)*(a + b + c - 2)*0.125; break;
                    case  7: M =  (a - 1)*(b + 1)*(c + 1)*(a - b - c + 2)*0.125; break;
                    case  8: M = -(a*a - 1)*(b - 1)*(c - 1)*0.250; break;
                    case  9: M =  (a + 1)*(b*b - 1)*(c - 1)*0.250; break;
                    case 10: M =  (a*a - 1)*(b + 1)*(c - 1)*0.250; break;
                    case 11: M = -(a - 1)*(b*b - 1)*(c - 1)*0.250; break;
                    case 12: M =  (a*a - 1)*(b - 1)*(c + 1)*0.250; break;
                    case 13: M = -(a + 1)*(b*b - 1)*(c + 1)*0.250; break;
                    case 14: M = -(a*a - 1)*(b + 1)*(c + 1)*0.250; break;
                    case 15: M =  (a - 1)*(b*b - 1)*(c + 1)*0.250; break;
                    case 16: M = -(a - 1)*(b - 1)*(c*c - 1)*0.250; break;
                    case 17: M =  (a + 1)*(b - 1)*(c*c - 1)*0.250; break;
                    case 18: M = -(a + 1)*(b + 1)*(c*c - 1)*0.250; break;
                    case 19: M =  (a - 1)*(b + 1)*(c*c - 1)*0.250; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                double derivative_Mi_c = 0.;
                switch (indx_i)
                {
                    case  0: derivative_Mi_a = (b - 1)*(c - 1)*(2*a + b + c + 1)*0.125;  derivative_Mi_b = (a - 1)*(c - 1)*(a + 2*b + c + 1)*0.125;  derivative_Mi_c = (a - 1)*(b - 1)*(a + b + 2*c + 1)*0.125;  break;
                    case  1: derivative_Mi_a = (b - 1)*(c - 1)*(2*a - b - c - 1)*0.125;  derivative_Mi_b = (a + 1)*(c - 1)*(a - 2*b - c - 1)*0.125;  derivative_Mi_c = (a + 1)*(b - 1)*(a - b - 2*c - 1)*0.125;  break;
                    case  2: derivative_Mi_a = (b + 1)*(c - 1)*(-2*a - b + c + 1)*0.125; derivative_Mi_b = (a + 1)*(c - 1)*(-a - 2*b + c + 1)*0.125; derivative_Mi_c = (a + 1)*(b + 1)*(-a - b + 2*c + 1)*0.125; break;
                    case  3: derivative_Mi_a = (b + 1)*(c - 1)*(-2*a + b - c - 1)*0.125; derivative_Mi_b = (a - 1)*(c - 1)*(-a + 2*b - c - 1)*0.125; derivative_Mi_c = (a - 1)*(b + 1)*(-a + b - 2*c - 1)*0.125; break;
                    case  4: derivative_Mi_a = (b - 1)*(c + 1)*(-2*a - b + c - 1)*0.125; derivative_Mi_b = (a - 1)*(c + 1)*(-a - 2*b + c - 1)*0.125; derivative_Mi_c = (a - 1)*(b - 1)*(-a - b + 2*c - 1)*0.125; break;
                    case  5: derivative_Mi_a = (b - 1)*(c + 1)*(-2*a + b - c + 1)*0.125; derivative_Mi_b = (a + 1)*(c + 1)*(-a + 2*b - c + 1)*0.125; derivative_Mi_c = (a + 1)*(b - 1)*(-a + b - 2*c + 1)*0.125; break;
                    case  6: derivative_Mi_a = (b + 1)*(c + 1)*(2*a + b + c - 1)*0.125;  derivative_Mi_b = (a + 1)*(c + 1)*(a + 2*b + c - 1)*0.125;  derivative_Mi_c = (a + 1)*(b + 1)*(a + b + 2*c - 1)*0.125;  break;
                    case  7: derivative_Mi_a = (b + 1)*(c + 1)*(2*a - b - c + 1)*0.125;  derivative_Mi_b = (a - 1)*(c + 1)*(a - 2*b - c + 1)*0.125;  derivative_Mi_c = (a - 1)*(b + 1)*(a - b - 2*c + 1)*0.125;  break;
                    case  8: derivative_Mi_a = -a*(b - 1)*(c - 1)*0.500;                 derivative_Mi_b = -(a*a - 1)*(c - 1)*0.250;                 derivative_Mi_c = -(a*a - 1)*(b - 1)*0.250;                 break;
                    case  9: derivative_Mi_a = (b*b - 1)*(c - 1)*0.250;                  derivative_Mi_b = b*(a + 1)*(c - 1)*0.500;                  derivative_Mi_c = (a + 1)*(b*b - 1)*0.250;                  break;
                    case 10: derivative_Mi_a = a*(b + 1)*(c - 1)*0.500;                  derivative_Mi_b = (a*a - 1)*(c - 1)*0.250;                  derivative_Mi_c = (a*a - 1)*(b + 1)*0.250;                  break;
                    case 11: derivative_Mi_a = -(b*b - 1)*(c - 1)*0.250;                 derivative_Mi_b = -b*(a - 1)*(c - 1)*0.500;                 derivative_Mi_c = -(a - 1)*(b*b - 1)*0.250;                 break;
                    case 12: derivative_Mi_a = a*(b - 1)*(c + 1)*0.500;                  derivative_Mi_b = (a*a - 1)*(c + 1)*0.250;                  derivative_Mi_c = (a*a - 1)*(b - 1)*0.250;                  break;
                    case 13: derivative_Mi_a = -(b*b - 1)*(c + 1)*0.250;                 derivative_Mi_b = -b*(a + 1)*(c + 1)*0.500;                 derivative_Mi_c = -(a + 1)*(b*b - 1)*0.250;                 break;
                    case 14: derivative_Mi_a = -a*(b + 1)*(c + 1)*0.500;                 derivative_Mi_b = -(a*a - 1)*(c + 1)*0.250;                 derivative_Mi_c = -(a*a - 1)*(b + 1)*0.250;                 break;
                    case 15: derivative_Mi_a = (b*b - 1)*(c + 1)*0.250;                  derivative_Mi_b = b*(a - 1)*(c + 1)*0.500;                  derivative_Mi_c = (a - 1)*(b*b - 1)*0.250;                  break;
                    case 16: derivative_Mi_a = -(b - 1)*(c*c - 1)*0.250;                 derivative_Mi_b = -(a - 1)*(c*c - 1)*0.250;                 derivative_Mi_c = -c*(a - 1)*(b - 1)*0.500;                 break;
                    case 17: derivative_Mi_a = (b - 1)*(c*c - 1)*0.250;                  derivative_Mi_b = (a + 1)*(c*c - 1)*0.250;                  derivative_Mi_c = c*(a + 1)*(b - 1)*0.500;                  break;
                    case 18: derivative_Mi_a = -(b + 1)*(c*c - 1)*0.250;                 derivative_Mi_b = -(a + 1)*(c*c - 1)*0.250;                 derivative_Mi_c = -c*(a + 1)*(b + 1)*0.500;                 break;
                    case 19: derivative_Mi_a = (b + 1)*(c*c - 1)*0.250;                  derivative_Mi_b = (a - 1)*(c*c - 1)*0.250;                  derivative_Mi_c = c*(a - 1)*(b + 1)*0.500;                  break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Mi_abc_vec;
                derivative_Mi_abc_vec << derivative_Mi_a, derivative_Mi_b, derivative_Mi_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Mi_xyz_vec = derivative_Mi_abc_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xyz_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xyz_vec.coeffRef(1);
                double derivative_Mi_z = derivative_Mi_xyz_vec.coeffRef(2);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);
                derivative_Mi_z_part_mli_vec.push_back(derivative_Mi_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);
            derivative_Mi_z_part_ml_vec.push_back(derivative_Mi_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_line_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_line_vec.push_back(N_part_ml_vec);
        Mi_line_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_line_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_line_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_line_vec.push_back(derivative_Ni_z_part_ml_vec);
        derivative_Mi_x_line_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_line_vec.push_back(derivative_Mi_y_part_ml_vec);
        derivative_Mi_z_line_vec.push_back(derivative_Mi_z_part_ml_vec);
 
    }

    // quadratic integration points

    // weights for integration
    weight_quad_vec = {
        0.1714677641, 0.1714677641, 0.1714677641, 0.1714677641, 0.1714677641, 0.1714677641, 0.1714677641, 0.1714677641,
        0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225,
        0.2743484225, 0.2743484225, 0.2743484225, 0.2743484225, 0.4389574760, 0.4389574760, 0.4389574760, 0.4389574760,
        0.4389574760, 0.4389574760, 0.7023319616
    };

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_quad_arr[27] = {
        -0.7745966692, +0.7745966692, -0.7745966692, +0.7745966692, -0.7745966692, +0.7745966692, -0.7745966692, +0.7745966692,
         0,            -0.7745966692, +0.7745966692, 0,             -0.7745966692, +0.7745966692, -0.7745966692, +0.7745966692,
         0,            -0.7745966692, +0.7745966692, 0,              0,             0,            -0.7745966692, +0.7745966692,
         0,             0,             0
    };
    double b_quad_arr[27] = {
        -0.7745966692, -0.7745966692, +0.7745966692, +0.7745966692, -0.7745966692, -0.7745966692, +0.7745966692, +0.7745966692,
        -0.7745966692,  0,             0,            +0.7745966692, -0.7745966692, -0.7745966692, +0.7745966692, +0.7745966692,
        -0.7745966692,  0,             0,            +0.7745966692,  0,            -0.7745966692,  0,             0,
        +0.7745966692,  0,             0
    };
    double c_quad_arr[27] = {
        -0.7745966692, -0.7745966692, -0.7745966692, -0.7745966692, +0.7745966692, +0.7745966692, +0.7745966692, +0.7745966692,
        -0.7745966692, -0.7745966692, -0.7745966692, -0.7745966692,  0,             0,             0,             0,
        +0.7745966692, +0.7745966692, +0.7745966692, +0.7745966692, -0.7745966692,  0,             0,             0,
         0,            +0.7745966692,  0
    };

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Ni_z_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;
        VectorDouble2D derivative_Mi_z_part_ml_vec;

        // get points around element
        int p00_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 0];
        int p01_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 1];
        int p02_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 2];
        int p03_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 3];
        int p04_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 4];
        int p05_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 5];
        int p06_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 6];
        int p07_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 7];
        int p08_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 8];
        int p09_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][ 9];
        int p10_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][10];
        int p11_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][11];
        int p12_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][12];
        int p13_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][13];
        int p14_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][14];
        int p15_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][15];
        int p16_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][16];
        int p17_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][17];
        int p18_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][18];
        int p19_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][19];
        int p00_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p00_pgid];
        int p01_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p01_pgid];
        int p02_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p02_pgid];
        int p03_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p03_pgid];
        int p04_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p04_pgid];
        int p05_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p05_pgid];
        int p06_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p06_pgid];
        int p07_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p07_pgid];
        int p08_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p08_pgid];
        int p09_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p09_pgid];
        int p10_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p10_pgid];
        int p11_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p11_pgid];
        int p12_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p12_pgid];
        int p13_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p13_pgid];
        int p14_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p14_pgid];
        int p15_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p15_pgid];
        int p16_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p16_pgid];
        int p17_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p17_pgid];
        int p18_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p18_pgid];
        int p19_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p19_pgid];

        // get x values of points
        double x00 = domain_quad_ptr->point_position_x_vec[p00_pdid];
        double x01 = domain_quad_ptr->point_position_x_vec[p01_pdid];
        double x02 = domain_quad_ptr->point_position_x_vec[p02_pdid];
        double x03 = domain_quad_ptr->point_position_x_vec[p03_pdid];
        double x04 = domain_quad_ptr->point_position_x_vec[p04_pdid];
        double x05 = domain_quad_ptr->point_position_x_vec[p05_pdid];
        double x06 = domain_quad_ptr->point_position_x_vec[p06_pdid];
        double x07 = domain_quad_ptr->point_position_x_vec[p07_pdid];
        double x08 = domain_quad_ptr->point_position_x_vec[p08_pdid];
        double x09 = domain_quad_ptr->point_position_x_vec[p09_pdid];
        double x10 = domain_quad_ptr->point_position_x_vec[p10_pdid];
        double x11 = domain_quad_ptr->point_position_x_vec[p11_pdid];
        double x12 = domain_quad_ptr->point_position_x_vec[p12_pdid];
        double x13 = domain_quad_ptr->point_position_x_vec[p13_pdid];
        double x14 = domain_quad_ptr->point_position_x_vec[p14_pdid];
        double x15 = domain_quad_ptr->point_position_x_vec[p15_pdid];
        double x16 = domain_quad_ptr->point_position_x_vec[p16_pdid];
        double x17 = domain_quad_ptr->point_position_x_vec[p17_pdid];
        double x18 = domain_quad_ptr->point_position_x_vec[p18_pdid];
        double x19 = domain_quad_ptr->point_position_x_vec[p19_pdid];

        // get y values of points
        double y00 = domain_quad_ptr->point_position_y_vec[p00_pdid];
        double y01 = domain_quad_ptr->point_position_y_vec[p01_pdid];
        double y02 = domain_quad_ptr->point_position_y_vec[p02_pdid];
        double y03 = domain_quad_ptr->point_position_y_vec[p03_pdid];
        double y04 = domain_quad_ptr->point_position_y_vec[p04_pdid];
        double y05 = domain_quad_ptr->point_position_y_vec[p05_pdid];
        double y06 = domain_quad_ptr->point_position_y_vec[p06_pdid];
        double y07 = domain_quad_ptr->point_position_y_vec[p07_pdid];
        double y08 = domain_quad_ptr->point_position_y_vec[p08_pdid];
        double y09 = domain_quad_ptr->point_position_y_vec[p09_pdid];
        double y10 = domain_quad_ptr->point_position_y_vec[p10_pdid];
        double y11 = domain_quad_ptr->point_position_y_vec[p11_pdid];
        double y12 = domain_quad_ptr->point_position_y_vec[p12_pdid];
        double y13 = domain_quad_ptr->point_position_y_vec[p13_pdid];
        double y14 = domain_quad_ptr->point_position_y_vec[p14_pdid];
        double y15 = domain_quad_ptr->point_position_y_vec[p15_pdid];
        double y16 = domain_quad_ptr->point_position_y_vec[p16_pdid];
        double y17 = domain_quad_ptr->point_position_y_vec[p17_pdid];
        double y18 = domain_quad_ptr->point_position_y_vec[p18_pdid];
        double y19 = domain_quad_ptr->point_position_y_vec[p19_pdid];

        // get z values of points
        double z00 = domain_quad_ptr->point_position_z_vec[p00_pdid];
        double z01 = domain_quad_ptr->point_position_z_vec[p01_pdid];
        double z02 = domain_quad_ptr->point_position_z_vec[p02_pdid];
        double z03 = domain_quad_ptr->point_position_z_vec[p03_pdid];
        double z04 = domain_quad_ptr->point_position_z_vec[p04_pdid];
        double z05 = domain_quad_ptr->point_position_z_vec[p05_pdid];
        double z06 = domain_quad_ptr->point_position_z_vec[p06_pdid];
        double z07 = domain_quad_ptr->point_position_z_vec[p07_pdid];
        double z08 = domain_quad_ptr->point_position_z_vec[p08_pdid];
        double z09 = domain_quad_ptr->point_position_z_vec[p09_pdid];
        double z10 = domain_quad_ptr->point_position_z_vec[p10_pdid];
        double z11 = domain_quad_ptr->point_position_z_vec[p11_pdid];
        double z12 = domain_quad_ptr->point_position_z_vec[p12_pdid];
        double z13 = domain_quad_ptr->point_position_z_vec[p13_pdid];
        double z14 = domain_quad_ptr->point_position_z_vec[p14_pdid];
        double z15 = domain_quad_ptr->point_position_z_vec[p15_pdid];
        double z16 = domain_quad_ptr->point_position_z_vec[p16_pdid];
        double z17 = domain_quad_ptr->point_position_z_vec[p17_pdid];
        double z18 = domain_quad_ptr->point_position_z_vec[p18_pdid];
        double z19 = domain_quad_ptr->point_position_z_vec[p19_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 27; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Mi_z_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;
            VectorDouble derivative_Ni_z_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_quad_arr[indx_l];
            double b = b_quad_arr[indx_l];
            double c = c_quad_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = -a*x08*(b - 1)*(c - 1)*0.500 + a*x10*(b + 1)*(c - 1)*0.500 + a*x12*(b - 1)*(c + 1)*0.500 - a*x14*(b + 1)*(c + 1)*0.500 + x00*(a - 1)*(b - 1)*(c - 1)*0.125 + x00*(b - 1)*(c - 1)*(a + b + c + 2)*0.125 + x01*(a + 1)*(b - 1)*(c - 1)*0.125 - x01*(b - 1)*(c - 1)*(-a + b + c + 2)*0.125 - x02*(a + 1)*(b + 1)*(c - 1)*0.125 - x02*(b + 1)*(c - 1)*(a + b - c - 2)*0.125 - x03*(a - 1)*(b + 1)*(c - 1)*0.125 - x03*(b + 1)*(c - 1)*(a - b + c + 2)*0.125 - x04*(a - 1)*(b - 1)*(c + 1)*0.125 - x04*(b - 1)*(c + 1)*(a + b - c + 2)*0.125 - x05*(a + 1)*(b - 1)*(c + 1)*0.125 - x05*(b - 1)*(c + 1)*(a - b + c - 2)*0.125 + x06*(a + 1)*(b + 1)*(c + 1)*0.125 + x06*(b + 1)*(c + 1)*(a + b + c - 2)*0.125 + x07*(a - 1)*(b + 1)*(c + 1)*0.125 + x07*(b + 1)*(c + 1)*(a - b - c + 2)*0.125 + x09*(b*b - 1)*(c - 1)*0.250 - x11*(b*b - 1)*(c - 1)*0.250 - x13*(b*b - 1)*(c + 1)*0.250 + x15*(b*b - 1)*(c + 1)*0.250 - x16*(b - 1)*(c*c - 1)*0.250 + x17*(b - 1)*(c*c - 1)*0.250 - x18*(b + 1)*(c*c - 1)*0.250 + x19*(b + 1)*(c*c - 1)*0.250;
            double derivative_x_b =  b*x09*(a + 1)*(c - 1)*0.500 - b*x11*(a - 1)*(c - 1)*0.500 - b*x13*(a + 1)*(c + 1)*0.500 + b*x15*(a - 1)*(c + 1)*0.500 + x00*(a - 1)*(b - 1)*(c - 1)*0.125 + x00*(a - 1)*(c - 1)*(a + b + c + 2)*0.125 - x01*(a + 1)*(b - 1)*(c - 1)*0.125 - x01*(a + 1)*(c - 1)*(-a + b + c + 2)*0.125 - x02*(a + 1)*(b + 1)*(c - 1)*0.125 - x02*(a + 1)*(c - 1)*(a + b - c - 2)*0.125 + x03*(a - 1)*(b + 1)*(c - 1)*0.125 - x03*(a - 1)*(c - 1)*(a - b + c + 2)*0.125 - x04*(a - 1)*(b - 1)*(c + 1)*0.125 - x04*(a - 1)*(c + 1)*(a + b - c + 2)*0.125 + x05*(a + 1)*(b - 1)*(c + 1)*0.125 - x05*(a + 1)*(c + 1)*(a - b + c - 2)*0.125 + x06*(a + 1)*(b + 1)*(c + 1)*0.125 + x06*(a + 1)*(c + 1)*(a + b + c - 2)*0.125 - x07*(a - 1)*(b + 1)*(c + 1)*0.125 + x07*(a - 1)*(c + 1)*(a - b - c + 2)*0.125 - x08*(a*a - 1)*(c - 1)*0.250 + x10*(a*a - 1)*(c - 1)*0.250 + x12*(a*a - 1)*(c + 1)*0.250 - x14*(a*a - 1)*(c + 1)*0.250 - x16*(a - 1)*(c*c - 1)*0.250 + x17*(a + 1)*(c*c - 1)*0.250 - x18*(a + 1)*(c*c - 1)*0.250 + x19*(a - 1)*(c*c - 1)*0.250;
            double derivative_x_c = -c*x16*(a - 1)*(b - 1)*0.500 + c*x17*(a + 1)*(b - 1)*0.500 - c*x18*(a + 1)*(b + 1)*0.500 + c*x19*(a - 1)*(b + 1)*0.500 + x00*(a - 1)*(b - 1)*(c - 1)*0.125 + x00*(a - 1)*(b - 1)*(a + b + c + 2)*0.125 - x01*(a + 1)*(b - 1)*(c - 1)*0.125 - x01*(a + 1)*(b - 1)*(-a + b + c + 2)*0.125 + x02*(a + 1)*(b + 1)*(c - 1)*0.125 - x02*(a + 1)*(b + 1)*(a + b - c - 2)*0.125 - x03*(a - 1)*(b + 1)*(c - 1)*0.125 - x03*(a - 1)*(b + 1)*(a - b + c + 2)*0.125 + x04*(a - 1)*(b - 1)*(c + 1)*0.125 - x04*(a - 1)*(b - 1)*(a + b - c + 2)*0.125 - x05*(a + 1)*(b - 1)*(c + 1)*0.125 - x05*(a + 1)*(b - 1)*(a - b + c - 2)*0.125 + x06*(a + 1)*(b + 1)*(c + 1)*0.125 + x06*(a + 1)*(b + 1)*(a + b + c - 2)*0.125 - x07*(a - 1)*(b + 1)*(c + 1)*0.125 + x07*(a - 1)*(b + 1)*(a - b - c + 2)*0.125 - x08*(a*a - 1)*(b - 1)*0.250 + x09*(a + 1)*(b*b - 1)*0.250 + x10*(a*a - 1)*(b + 1)*0.250 - x11*(a - 1)*(b*b - 1)*0.250 + x12*(a*a - 1)*(b - 1)*0.250 - x13*(a + 1)*(b*b - 1)*0.250 - x14*(a*a - 1)*(b + 1)*0.250 + x15*(a - 1)*(b*b - 1)*0.250;
            double derivative_y_a = -a*y08*(b - 1)*(c - 1)*0.500 + a*y10*(b + 1)*(c - 1)*0.500 + a*y12*(b - 1)*(c + 1)*0.500 - a*y14*(b + 1)*(c + 1)*0.500 + y00*(a - 1)*(b - 1)*(c - 1)*0.125 + y00*(b - 1)*(c - 1)*(a + b + c + 2)*0.125 + y01*(a + 1)*(b - 1)*(c - 1)*0.125 - y01*(b - 1)*(c - 1)*(-a + b + c + 2)*0.125 - y02*(a + 1)*(b + 1)*(c - 1)*0.125 - y02*(b + 1)*(c - 1)*(a + b - c - 2)*0.125 - y03*(a - 1)*(b + 1)*(c - 1)*0.125 - y03*(b + 1)*(c - 1)*(a - b + c + 2)*0.125 - y04*(a - 1)*(b - 1)*(c + 1)*0.125 - y04*(b - 1)*(c + 1)*(a + b - c + 2)*0.125 - y05*(a + 1)*(b - 1)*(c + 1)*0.125 - y05*(b - 1)*(c + 1)*(a - b + c - 2)*0.125 + y06*(a + 1)*(b + 1)*(c + 1)*0.125 + y06*(b + 1)*(c + 1)*(a + b + c - 2)*0.125 + y07*(a - 1)*(b + 1)*(c + 1)*0.125 + y07*(b + 1)*(c + 1)*(a - b - c + 2)*0.125 + y09*(b*b - 1)*(c - 1)*0.250 - y11*(b*b - 1)*(c - 1)*0.250 - y13*(b*b - 1)*(c + 1)*0.250 + y15*(b*b - 1)*(c + 1)*0.250 - y16*(b - 1)*(c*c - 1)*0.250 + y17*(b - 1)*(c*c - 1)*0.250 - y18*(b + 1)*(c*c - 1)*0.250 + y19*(b + 1)*(c*c - 1)*0.250;
            double derivative_y_b =  b*y09*(a + 1)*(c - 1)*0.500 - b*y11*(a - 1)*(c - 1)*0.500 - b*y13*(a + 1)*(c + 1)*0.500 + b*y15*(a - 1)*(c + 1)*0.500 + y00*(a - 1)*(b - 1)*(c - 1)*0.125 + y00*(a - 1)*(c - 1)*(a + b + c + 2)*0.125 - y01*(a + 1)*(b - 1)*(c - 1)*0.125 - y01*(a + 1)*(c - 1)*(-a + b + c + 2)*0.125 - y02*(a + 1)*(b + 1)*(c - 1)*0.125 - y02*(a + 1)*(c - 1)*(a + b - c - 2)*0.125 + y03*(a - 1)*(b + 1)*(c - 1)*0.125 - y03*(a - 1)*(c - 1)*(a - b + c + 2)*0.125 - y04*(a - 1)*(b - 1)*(c + 1)*0.125 - y04*(a - 1)*(c + 1)*(a + b - c + 2)*0.125 + y05*(a + 1)*(b - 1)*(c + 1)*0.125 - y05*(a + 1)*(c + 1)*(a - b + c - 2)*0.125 + y06*(a + 1)*(b + 1)*(c + 1)*0.125 + y06*(a + 1)*(c + 1)*(a + b + c - 2)*0.125 - y07*(a - 1)*(b + 1)*(c + 1)*0.125 + y07*(a - 1)*(c + 1)*(a - b - c + 2)*0.125 - y08*(a*a - 1)*(c - 1)*0.250 + y10*(a*a - 1)*(c - 1)*0.250 + y12*(a*a - 1)*(c + 1)*0.250 - y14*(a*a - 1)*(c + 1)*0.250 - y16*(a - 1)*(c*c - 1)*0.250 + y17*(a + 1)*(c*c - 1)*0.250 - y18*(a + 1)*(c*c - 1)*0.250 + y19*(a - 1)*(c*c - 1)*0.250;
            double derivative_y_c = -c*y16*(a - 1)*(b - 1)*0.500 + c*y17*(a + 1)*(b - 1)*0.500 - c*y18*(a + 1)*(b + 1)*0.500 + c*y19*(a - 1)*(b + 1)*0.500 + y00*(a - 1)*(b - 1)*(c - 1)*0.125 + y00*(a - 1)*(b - 1)*(a + b + c + 2)*0.125 - y01*(a + 1)*(b - 1)*(c - 1)*0.125 - y01*(a + 1)*(b - 1)*(-a + b + c + 2)*0.125 + y02*(a + 1)*(b + 1)*(c - 1)*0.125 - y02*(a + 1)*(b + 1)*(a + b - c - 2)*0.125 - y03*(a - 1)*(b + 1)*(c - 1)*0.125 - y03*(a - 1)*(b + 1)*(a - b + c + 2)*0.125 + y04*(a - 1)*(b - 1)*(c + 1)*0.125 - y04*(a - 1)*(b - 1)*(a + b - c + 2)*0.125 - y05*(a + 1)*(b - 1)*(c + 1)*0.125 - y05*(a + 1)*(b - 1)*(a - b + c - 2)*0.125 + y06*(a + 1)*(b + 1)*(c + 1)*0.125 + y06*(a + 1)*(b + 1)*(a + b + c - 2)*0.125 - y07*(a - 1)*(b + 1)*(c + 1)*0.125 + y07*(a - 1)*(b + 1)*(a - b - c + 2)*0.125 - y08*(a*a - 1)*(b - 1)*0.250 + y09*(a + 1)*(b*b - 1)*0.250 + y10*(a*a - 1)*(b + 1)*0.250 - y11*(a - 1)*(b*b - 1)*0.250 + y12*(a*a - 1)*(b - 1)*0.250 - y13*(a + 1)*(b*b - 1)*0.250 - y14*(a*a - 1)*(b + 1)*0.250 + y15*(a - 1)*(b*b - 1)*0.250;
            double derivative_z_a = -a*z08*(b - 1)*(c - 1)*0.500 + a*z10*(b + 1)*(c - 1)*0.500 + a*z12*(b - 1)*(c + 1)*0.500 - a*z14*(b + 1)*(c + 1)*0.500 + z00*(a - 1)*(b - 1)*(c - 1)*0.125 + z00*(b - 1)*(c - 1)*(a + b + c + 2)*0.125 + z01*(a + 1)*(b - 1)*(c - 1)*0.125 - z01*(b - 1)*(c - 1)*(-a + b + c + 2)*0.125 - z02*(a + 1)*(b + 1)*(c - 1)*0.125 - z02*(b + 1)*(c - 1)*(a + b - c - 2)*0.125 - z03*(a - 1)*(b + 1)*(c - 1)*0.125 - z03*(b + 1)*(c - 1)*(a - b + c + 2)*0.125 - z04*(a - 1)*(b - 1)*(c + 1)*0.125 - z04*(b - 1)*(c + 1)*(a + b - c + 2)*0.125 - z05*(a + 1)*(b - 1)*(c + 1)*0.125 - z05*(b - 1)*(c + 1)*(a - b + c - 2)*0.125 + z06*(a + 1)*(b + 1)*(c + 1)*0.125 + z06*(b + 1)*(c + 1)*(a + b + c - 2)*0.125 + z07*(a - 1)*(b + 1)*(c + 1)*0.125 + z07*(b + 1)*(c + 1)*(a - b - c + 2)*0.125 + z09*(b*b - 1)*(c - 1)*0.250 - z11*(b*b - 1)*(c - 1)*0.250 - z13*(b*b - 1)*(c + 1)*0.250 + z15*(b*b - 1)*(c + 1)*0.250 - z16*(b - 1)*(c*c - 1)*0.250 + z17*(b - 1)*(c*c - 1)*0.250 - z18*(b + 1)*(c*c - 1)*0.250 + z19*(b + 1)*(c*c - 1)*0.250;
            double derivative_z_b =  b*z09*(a + 1)*(c - 1)*0.500 - b*z11*(a - 1)*(c - 1)*0.500 - b*z13*(a + 1)*(c + 1)*0.500 + b*z15*(a - 1)*(c + 1)*0.500 + z00*(a - 1)*(b - 1)*(c - 1)*0.125 + z00*(a - 1)*(c - 1)*(a + b + c + 2)*0.125 - z01*(a + 1)*(b - 1)*(c - 1)*0.125 - z01*(a + 1)*(c - 1)*(-a + b + c + 2)*0.125 - z02*(a + 1)*(b + 1)*(c - 1)*0.125 - z02*(a + 1)*(c - 1)*(a + b - c - 2)*0.125 + z03*(a - 1)*(b + 1)*(c - 1)*0.125 - z03*(a - 1)*(c - 1)*(a - b + c + 2)*0.125 - z04*(a - 1)*(b - 1)*(c + 1)*0.125 - z04*(a - 1)*(c + 1)*(a + b - c + 2)*0.125 + z05*(a + 1)*(b - 1)*(c + 1)*0.125 - z05*(a + 1)*(c + 1)*(a - b + c - 2)*0.125 + z06*(a + 1)*(b + 1)*(c + 1)*0.125 + z06*(a + 1)*(c + 1)*(a + b + c - 2)*0.125 - z07*(a - 1)*(b + 1)*(c + 1)*0.125 + z07*(a - 1)*(c + 1)*(a - b - c + 2)*0.125 - z08*(a*a - 1)*(c - 1)*0.250 + z10*(a*a - 1)*(c - 1)*0.250 + z12*(a*a - 1)*(c + 1)*0.250 - z14*(a*a - 1)*(c + 1)*0.250 - z16*(a - 1)*(c*c - 1)*0.250 + z17*(a + 1)*(c*c - 1)*0.250 - z18*(a + 1)*(c*c - 1)*0.250 + z19*(a - 1)*(c*c - 1)*0.250;
            double derivative_z_c = -c*z16*(a - 1)*(b - 1)*0.500 + c*z17*(a + 1)*(b - 1)*0.500 - c*z18*(a + 1)*(b + 1)*0.500 + c*z19*(a - 1)*(b + 1)*0.500 + z00*(a - 1)*(b - 1)*(c - 1)*0.125 + z00*(a - 1)*(b - 1)*(a + b + c + 2)*0.125 - z01*(a + 1)*(b - 1)*(c - 1)*0.125 - z01*(a + 1)*(b - 1)*(-a + b + c + 2)*0.125 + z02*(a + 1)*(b + 1)*(c - 1)*0.125 - z02*(a + 1)*(b + 1)*(a + b - c - 2)*0.125 - z03*(a - 1)*(b + 1)*(c - 1)*0.125 - z03*(a - 1)*(b + 1)*(a - b + c + 2)*0.125 + z04*(a - 1)*(b - 1)*(c + 1)*0.125 - z04*(a - 1)*(b - 1)*(a + b - c + 2)*0.125 - z05*(a + 1)*(b - 1)*(c + 1)*0.125 - z05*(a + 1)*(b - 1)*(a - b + c - 2)*0.125 + z06*(a + 1)*(b + 1)*(c + 1)*0.125 + z06*(a + 1)*(b + 1)*(a + b + c - 2)*0.125 - z07*(a - 1)*(b + 1)*(c + 1)*0.125 + z07*(a - 1)*(b + 1)*(a - b - c + 2)*0.125 - z08*(a*a - 1)*(b - 1)*0.250 + z09*(a + 1)*(b*b - 1)*0.250 + z10*(a*a - 1)*(b + 1)*0.250 - z11*(a - 1)*(b*b - 1)*0.250 + z12*(a*a - 1)*(b - 1)*0.250 - z13*(a + 1)*(b*b - 1)*0.250 - z14*(a*a - 1)*(b + 1)*0.250 + z15*(a - 1)*(b*b - 1)*0.250;

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

            // iterate for each test function
            for (int indx_i = 0; indx_i < 20; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case  0: M =  (a - 1)*(b - 1)*(c - 1)*(a + b + c + 2)*0.125; break;
                    case  1: M = -(a + 1)*(b - 1)*(c - 1)*(-a + b + c + 2)*0.125; break;
                    case  2: M = -(a + 1)*(b + 1)*(c - 1)*(a + b - c - 2)*0.125; break;
                    case  3: M = -(a - 1)*(b + 1)*(c - 1)*(a - b + c + 2)*0.125; break;
                    case  4: M = -(a - 1)*(b - 1)*(c + 1)*(a + b - c + 2)*0.125; break;
                    case  5: M = -(a + 1)*(b - 1)*(c + 1)*(a - b + c - 2)*0.125; break;
                    case  6: M =  (a + 1)*(b + 1)*(c + 1)*(a + b + c - 2)*0.125; break;
                    case  7: M =  (a - 1)*(b + 1)*(c + 1)*(a - b - c + 2)*0.125; break;
                    case  8: M = -(a*a - 1)*(b - 1)*(c - 1)*0.250; break;
                    case  9: M =  (a + 1)*(b*b - 1)*(c - 1)*0.250; break;
                    case 10: M =  (a*a - 1)*(b + 1)*(c - 1)*0.250; break;
                    case 11: M = -(a - 1)*(b*b - 1)*(c - 1)*0.250; break;
                    case 12: M =  (a*a - 1)*(b - 1)*(c + 1)*0.250; break;
                    case 13: M = -(a + 1)*(b*b - 1)*(c + 1)*0.250; break;
                    case 14: M = -(a*a - 1)*(b + 1)*(c + 1)*0.250; break;
                    case 15: M =  (a - 1)*(b*b - 1)*(c + 1)*0.250; break;
                    case 16: M = -(a - 1)*(b - 1)*(c*c - 1)*0.250; break;
                    case 17: M =  (a + 1)*(b - 1)*(c*c - 1)*0.250; break;
                    case 18: M = -(a + 1)*(b + 1)*(c*c - 1)*0.250; break;
                    case 19: M =  (a - 1)*(b + 1)*(c*c - 1)*0.250; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                double derivative_Mi_c = 0.;
                switch (indx_i)
                {
                    case  0: derivative_Mi_a = (b - 1)*(c - 1)*(2*a + b + c + 1)*0.125;  derivative_Mi_b = (a - 1)*(c - 1)*(a + 2*b + c + 1)*0.125;  derivative_Mi_c = (a - 1)*(b - 1)*(a + b + 2*c + 1)*0.125;  break;
                    case  1: derivative_Mi_a = (b - 1)*(c - 1)*(2*a - b - c - 1)*0.125;  derivative_Mi_b = (a + 1)*(c - 1)*(a - 2*b - c - 1)*0.125;  derivative_Mi_c = (a + 1)*(b - 1)*(a - b - 2*c - 1)*0.125;  break;
                    case  2: derivative_Mi_a = (b + 1)*(c - 1)*(-2*a - b + c + 1)*0.125; derivative_Mi_b = (a + 1)*(c - 1)*(-a - 2*b + c + 1)*0.125; derivative_Mi_c = (a + 1)*(b + 1)*(-a - b + 2*c + 1)*0.125; break;
                    case  3: derivative_Mi_a = (b + 1)*(c - 1)*(-2*a + b - c - 1)*0.125; derivative_Mi_b = (a - 1)*(c - 1)*(-a + 2*b - c - 1)*0.125; derivative_Mi_c = (a - 1)*(b + 1)*(-a + b - 2*c - 1)*0.125; break;
                    case  4: derivative_Mi_a = (b - 1)*(c + 1)*(-2*a - b + c - 1)*0.125; derivative_Mi_b = (a - 1)*(c + 1)*(-a - 2*b + c - 1)*0.125; derivative_Mi_c = (a - 1)*(b - 1)*(-a - b + 2*c - 1)*0.125; break;
                    case  5: derivative_Mi_a = (b - 1)*(c + 1)*(-2*a + b - c + 1)*0.125; derivative_Mi_b = (a + 1)*(c + 1)*(-a + 2*b - c + 1)*0.125; derivative_Mi_c = (a + 1)*(b - 1)*(-a + b - 2*c + 1)*0.125; break;
                    case  6: derivative_Mi_a = (b + 1)*(c + 1)*(2*a + b + c - 1)*0.125;  derivative_Mi_b = (a + 1)*(c + 1)*(a + 2*b + c - 1)*0.125;  derivative_Mi_c = (a + 1)*(b + 1)*(a + b + 2*c - 1)*0.125;  break;
                    case  7: derivative_Mi_a = (b + 1)*(c + 1)*(2*a - b - c + 1)*0.125;  derivative_Mi_b = (a - 1)*(c + 1)*(a - 2*b - c + 1)*0.125;  derivative_Mi_c = (a - 1)*(b + 1)*(a - b - 2*c + 1)*0.125;  break;
                    case  8: derivative_Mi_a = -a*(b - 1)*(c - 1)*0.500;                 derivative_Mi_b = -(a*a - 1)*(c - 1)*0.250;                 derivative_Mi_c = -(a*a - 1)*(b - 1)*0.250;                 break;
                    case  9: derivative_Mi_a = (b*b - 1)*(c - 1)*0.250;                  derivative_Mi_b = b*(a + 1)*(c - 1)*0.500;                  derivative_Mi_c = (a + 1)*(b*b - 1)*0.250;                  break;
                    case 10: derivative_Mi_a = a*(b + 1)*(c - 1)*0.500;                  derivative_Mi_b = (a*a - 1)*(c - 1)*0.250;                  derivative_Mi_c = (a*a - 1)*(b + 1)*0.250;                  break;
                    case 11: derivative_Mi_a = -(b*b - 1)*(c - 1)*0.250;                 derivative_Mi_b = -b*(a - 1)*(c - 1)*0.500;                 derivative_Mi_c = -(a - 1)*(b*b - 1)*0.250;                 break;
                    case 12: derivative_Mi_a = a*(b - 1)*(c + 1)*0.500;                  derivative_Mi_b = (a*a - 1)*(c + 1)*0.250;                  derivative_Mi_c = (a*a - 1)*(b - 1)*0.250;                  break;
                    case 13: derivative_Mi_a = -(b*b - 1)*(c + 1)*0.250;                 derivative_Mi_b = -b*(a + 1)*(c + 1)*0.500;                 derivative_Mi_c = -(a + 1)*(b*b - 1)*0.250;                 break;
                    case 14: derivative_Mi_a = -a*(b + 1)*(c + 1)*0.500;                 derivative_Mi_b = -(a*a - 1)*(c + 1)*0.250;                 derivative_Mi_c = -(a*a - 1)*(b + 1)*0.250;                 break;
                    case 15: derivative_Mi_a = (b*b - 1)*(c + 1)*0.250;                  derivative_Mi_b = b*(a - 1)*(c + 1)*0.500;                  derivative_Mi_c = (a - 1)*(b*b - 1)*0.250;                  break;
                    case 16: derivative_Mi_a = -(b - 1)*(c*c - 1)*0.250;                 derivative_Mi_b = -(a - 1)*(c*c - 1)*0.250;                 derivative_Mi_c = -c*(a - 1)*(b - 1)*0.500;                 break;
                    case 17: derivative_Mi_a = (b - 1)*(c*c - 1)*0.250;                  derivative_Mi_b = (a + 1)*(c*c - 1)*0.250;                  derivative_Mi_c = c*(a + 1)*(b - 1)*0.500;                  break;
                    case 18: derivative_Mi_a = -(b + 1)*(c*c - 1)*0.250;                 derivative_Mi_b = -(a + 1)*(c*c - 1)*0.250;                 derivative_Mi_c = -c*(a + 1)*(b + 1)*0.500;                 break;
                    case 19: derivative_Mi_a = (b + 1)*(c*c - 1)*0.250;                  derivative_Mi_b = (a - 1)*(c*c - 1)*0.250;                  derivative_Mi_c = c*(a - 1)*(b + 1)*0.500;                  break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector3d derivative_Mi_abc_vec;
                derivative_Mi_abc_vec << derivative_Mi_a, derivative_Mi_b, derivative_Mi_c;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector3d derivative_Mi_xyz_vec = derivative_Mi_abc_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xyz_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xyz_vec.coeffRef(1);
                double derivative_Mi_z = derivative_Mi_xyz_vec.coeffRef(2);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);
                derivative_Mi_z_part_mli_vec.push_back(derivative_Mi_z);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Ni_z_part_ml_vec.push_back(derivative_Ni_z_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);
            derivative_Mi_z_part_ml_vec.push_back(derivative_Mi_z_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_quad_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_quad_vec.push_back(N_part_ml_vec);
        Mi_quad_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_quad_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_quad_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Ni_z_quad_vec.push_back(derivative_Ni_z_part_ml_vec);
        derivative_Mi_x_quad_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_quad_vec.push_back(derivative_Mi_y_part_ml_vec);
        derivative_Mi_z_quad_vec.push_back(derivative_Mi_z_part_ml_vec);
 
    }

}

}

#endif
