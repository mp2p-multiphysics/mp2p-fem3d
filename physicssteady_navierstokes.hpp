#ifndef PHYSICSSTEADY_NAVIERSTOKES
#define PHYSICSSTEADY_NAVIERSTOKES
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_0d.hpp"
#include "domain_2d.hpp"
#include "domain_3d.hpp"
#include "integral_taylorhood_3d.hpp"
#include "physicssteady_base.hpp"
#include "scalar_0d.hpp"
#include "scalar_2d.hpp"
#include "scalar_3d.hpp"
#include "variable_3d.hpp"
#include "variable_group.hpp"

namespace FEM3D
{

class PhysicsSteadyNavierStokes : public PhysicsSteadyBase
{
    /*

    Steady-state Navier-Stokes equation (NSE) and continuity equation.
    
    NSE x-component: rho * (u * du/dx + v * du/dy + w * du/dz) = -dp/dx + mu * (d^2(u)/dx^2 + d^2(u)/dy^2 + d^2(u)/dz^2) + rho * g_x
    NSE y-component: rho * (u * dv/dx + v * dv/dy + w * dv/dz) = -dp/dy + mu * (d^2(v)/dx^2 + d^2(v)/dy^2 + d^2(v)/dz^2) + rho * g_y
    NSE z-component: rho * (u * dw/dx + v * dw/dy + w * dw/dz) = -dp/dz + mu * (d^2(w)/dx^2 + d^2(w)/dy^2 + d^2(w)/dz^2) + rho * g_z
    Continuity: du/dx + dv/dy + dw/dz = 0

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_variablegroup : void
        Set variables used in this physics.
    set_domain : void
        Set scalars applied to 2D domains.
    set_boundary_velocity : void
        Set a velocity boundary condition along a 1D domain.
    set_boundary_pressure : void
        Set a pressure boundary condition along a 1D domain.

    */

    public:

    // variables
    VariableGroup* velocity_x_ptr;
    VariableGroup* velocity_y_ptr;
    VariableGroup* velocity_z_ptr;
    VariableGroup* pressure_ptr;

    // domain objects
    std::vector<Domain3D*> domain_line_ptr_vec;
    std::vector<Domain3D*> domain_quad_ptr_vec;
    std::vector<IntegralTaylorHood3D*> integral_ptr_vec;
    std::vector<Scalar3D*> density_ptr_vec;
    std::vector<Scalar3D*> viscosity_ptr_vec;
    std::vector<Scalar3D*> force_x_ptr_vec;
    std::vector<Scalar3D*> force_y_ptr_vec;
    std::vector<Scalar3D*> force_z_ptr_vec;

    // boundary objects - velocity
    std::vector<Domain2D*> velocity_domain_ptr_vec;
    std::vector<Scalar2D*> velocity_x_ptr_vec;
    std::vector<Scalar2D*> velocity_y_ptr_vec;
    std::vector<Scalar2D*> velocity_z_ptr_vec;

    // boundary objects - pressure
    std::vector<Domain2D*> pressure_domain_ptr_vec;
    std::vector<Scalar2D*> pressure_ptr_vec;

    // boundary objects - pressure point
    std::vector<Domain0D*> pressure_point_domain_ptr_vec;
    std::vector<Scalar0D*> pressure_point_ptr_vec;

    // vectors of objects to update
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<Scalar3D*> scalar3d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec);
    void set_variablegroup(VariableGroup &velocity_x_in, VariableGroup &velocity_y_in, VariableGroup &velocity_z_in, VariableGroup &pressure_in);
    void set_domain(Domain3D &domain_line_in, Domain3D &domain_quad_in, IntegralTaylorHood3D &integral_in, Scalar3D &density_in, Scalar3D &viscosity_in, Scalar3D &force_x_in, Scalar3D &force_y_in, Scalar3D &force_z_in);
    void set_boundary_velocity(Domain2D &domain_quad_in, Scalar2D &velocity_x_in, Scalar2D &velocity_y_in, Scalar2D &velocity_z_in);
    void set_boundary_pressure(Domain2D &domain_quad_in, Scalar2D &pressure_in); 
    void set_boundary_pressure(Domain0D &domain_in, Scalar0D &pressure_in); 

    // getter and setter functions
    void set_start_row(int start_row_in) {start_row = start_row_in;}
    int get_start_row() {return start_row;}
    std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    std::vector<Scalar3D*> get_scalar3d_ptr_vec() {return scalar3d_ptr_vec;}
    std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsSteadyNavierStokes() {}

    private:

    void matrix_fill_domain
    (
        std::vector<EigenTriplet> &delta_a_triplet_vec, EigenVector &b_vec, EigenVector &x_vec,
        Domain3D *domain_line_ptr, Domain3D *domain_quad_ptr, IntegralTaylorHood3D *integral_ptr,
        Scalar3D *density_ptr, Scalar3D *viscosity_ptr, Scalar3D *force_x_ptr, Scalar3D *force_y_ptr, Scalar3D *force_z_ptr
    );
    void matrix_fill_velocity_clear
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr
    );
    void matrix_fill_pressure_clear
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr
    );
    void matrix_fill_pressure_point_clear
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain0D *domain_ptr
    );
    void matrix_fill_velocity
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr,
        Scalar2D *velocity_x_ptr, Scalar2D *velocity_y_ptr, Scalar2D *velocity_z_ptr
    );
    void matrix_fill_pressure
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr,
        Scalar2D *pressure_ptr
    );
    void matrix_fill_pressure_point
    (
        EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
        Domain0D *domain_ptr,
        Scalar0D *pressure_ptr
    );

};

void PhysicsSteadyNavierStokes::set_variablegroup(VariableGroup &velocity_x_in, VariableGroup &velocity_y_in, VariableGroup &velocity_z_in, VariableGroup &pressure_in)
{
    /*
    
    Set variables used in this physics.

    Arguments
    =========
    velocity_x_in : VariableGroup
        u in the NSE and continuity equation.
        Must be defined over a quadratic domain.
    velocity_y_in : VariableGroup
        v in the NSE and continuity equation.
        Must be defined over a quadratic domain.
    velocity_z_in : VariableGroup
        w in the NSE and continuity equation.
        Must be defined over a quadratic domain.
    pressure_in : VariableGroup
        p in the NSE.
        Must be defined over a linear domain.

    Returns
    =======
    (none)
    
    */

    // set variable groups
    velocity_x_ptr = &velocity_x_in;
    velocity_y_ptr = &velocity_y_in;
    velocity_z_ptr = &velocity_z_in;
    pressure_ptr = &pressure_in;

    // add to vector of variable groups
    variablegroup_ptr_vec.push_back(&velocity_x_in);
    variablegroup_ptr_vec.push_back(&velocity_y_in);
    variablegroup_ptr_vec.push_back(&velocity_z_in);
    variablegroup_ptr_vec.push_back(&pressure_in);

}

void PhysicsSteadyNavierStokes::set_domain(Domain3D &domain_line_in, Domain3D &domain_quad_in, IntegralTaylorHood3D &integral_in, Scalar3D &density_in, Scalar3D &viscosity_in, Scalar3D &force_x_in, Scalar3D &force_y_in, Scalar3D &force_z_in)
{
    /*
    
    Set scalars applied to 2D domains.

    Arguments
    =========
    domain_line_in : Domain3D
        Linear domain that this physics applies to.
    domain_quad_in : Domain3D
        Quadratic domain that this physics applies to.
    integral_in : IntegralTaylorHood3D
        Test function integrals over the domains.
    density_in : Scalar3D
        rho in NSE.
    viscosity_in : Scalar3D
        mu in NSE.
    force_x_in : Scalar3D
        g_x in NSE.
    force_y_in : Scalar3D
        g_y in NSE.
    force_z_in : Scalar3D
        g_z in NSE.

    Returns
    =======
    (none)

    Notes
    =====
    Coinciding points in the linear and quadratic domains must have the same IDs.

    */

    // add to vector of domain objects
    domain_line_ptr_vec.push_back(&domain_line_in);
    domain_quad_ptr_vec.push_back(&domain_quad_in);
    integral_ptr_vec.push_back(&integral_in);
    density_ptr_vec.push_back(&density_in);
    viscosity_ptr_vec.push_back(&viscosity_in);
    force_x_ptr_vec.push_back(&force_x_in);
    force_y_ptr_vec.push_back(&force_y_in);
    force_z_ptr_vec.push_back(&force_z_in);

    // add to vector of scalar3d objects
    scalar3d_ptr_vec.push_back(&density_in);
    scalar3d_ptr_vec.push_back(&viscosity_in);
    scalar3d_ptr_vec.push_back(&force_x_in);
    scalar3d_ptr_vec.push_back(&force_y_in);
    scalar3d_ptr_vec.push_back(&force_z_in);

    // calculate linear domain integrals
    integral_in.evaluate_integral_Ni_derivative_Mj_x_line();
    integral_in.evaluate_integral_Ni_derivative_Mj_y_line();
    integral_in.evaluate_integral_Ni_derivative_Mj_z_line();
    integral_in.evaluate_integral_Ni_Nj_line();

    // calculate quadratic domain integrals
    integral_in.evaluate_integral_Mi_quad();
    integral_in.evaluate_integral_Mi_derivative_Nj_x_quad();
    integral_in.evaluate_integral_Mi_derivative_Nj_y_quad();
    integral_in.evaluate_integral_Mi_derivative_Nj_z_quad();
    integral_in.evaluate_integral_div_Mi_dot_div_Mj_quad();
    integral_in.evaluate_integral_Mi_Mk_derivative_Mj_x_quad();
    integral_in.evaluate_integral_Mi_Mk_derivative_Mj_y_quad();
    integral_in.evaluate_integral_Mi_Mk_derivative_Mj_z_quad();

}

void PhysicsSteadyNavierStokes::set_boundary_velocity(Domain2D &domain_quad_in, Scalar2D &velocity_x_in, Scalar2D &velocity_y_in, Scalar2D &velocity_z_in)
{
    /*
    
    Set a velocity boundary condition along a 1D domain.

    Arguments
    =========
    domain_quad_in : Domain2D
        Quadratic domain that this boundary condition applies to.
    velocity_x_in : Scalar2D
        x velocity prescribed by the boundary condition.
    velocity_y_in : Scalar2D
        y velocity prescribed by the boundary condition.
    velocity_z_in : Scalar2D
        z velocity prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    velocity_domain_ptr_vec.push_back(&domain_quad_in);
    velocity_x_ptr_vec.push_back(&velocity_x_in);
    velocity_y_ptr_vec.push_back(&velocity_y_in);
    velocity_z_ptr_vec.push_back(&velocity_z_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&velocity_x_in);
    scalar2d_ptr_vec.push_back(&velocity_y_in);
    scalar2d_ptr_vec.push_back(&velocity_z_in);

}

void PhysicsSteadyNavierStokes::set_boundary_pressure(Domain2D &domain_line_in, Scalar2D &pressure_in)
{
    /*
    
    Set a pressure boundary condition along a 1D domain.

    Arguments
    =========
    domain_line_in : Domain2D
        Linear domain that this boundary condition applies to.
    pressure_in : Scalar2D
        Pressure prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    pressure_domain_ptr_vec.push_back(&domain_line_in);
    pressure_ptr_vec.push_back(&pressure_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&pressure_in);

}

void PhysicsSteadyNavierStokes::set_boundary_pressure(Domain0D &domain_in, Scalar0D &pressure_in)
{
    /*
    
    Set a pressure boundary condition along a 0D domain.

    Arguments
    =========
    domain_in : Domain2D
        Point domain that this boundary condition applies to.
    pressure_in : Scalar2D
        Pressure prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    pressure_point_domain_ptr_vec.push_back(&domain_in);
    pressure_point_ptr_vec.push_back(&pressure_in);

}

void PhysicsSteadyNavierStokes::matrix_fill
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec
)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : EigenSparseMatrix
        A in Ax = b.
    b_vec : EigenVector
        b in Ax = b.
    x_vec : EigenVector
        x in Ax = b.

    Returns
    =======
    (none)

    */

    // represent matrix as triplets for performance
    std::vector<EigenTriplet> delta_a_triplet_vec;
    int num_equation = a_mat.rows();
    delta_a_triplet_vec.reserve(20*num_equation); // estimated number of entries

   // iterate through each domain
   for (int indx_d = 0; indx_d < domain_line_ptr_vec.size(); indx_d++)
   {

        // subset domain objects
        Domain3D *domain_line_ptr = domain_line_ptr_vec[indx_d];
        Domain3D *domain_quad_ptr = domain_quad_ptr_vec[indx_d];
        IntegralTaylorHood3D *integral_ptr = integral_ptr_vec[indx_d];
        Scalar3D *density_ptr = density_ptr_vec[indx_d];
        Scalar3D *viscosity_ptr = viscosity_ptr_vec[indx_d];
        Scalar3D *force_x_ptr = force_x_ptr_vec[indx_d];
        Scalar3D *force_y_ptr = force_y_ptr_vec[indx_d];
        Scalar3D *force_z_ptr = force_z_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(delta_a_triplet_vec, b_vec, x_vec, domain_line_ptr, domain_quad_ptr, integral_ptr, density_ptr, viscosity_ptr, force_x_ptr, force_y_ptr, force_z_ptr);

   }

    // convert triplet vector to sparse matrix
    EigenSparseMatrix delta_a_mat(num_equation, num_equation);
    delta_a_mat.setFromTriplets(delta_a_triplet_vec.begin(), delta_a_triplet_vec.end());
    a_mat += delta_a_mat;

    // clear equations with dirichlet boundary conditions
    for (int indx_d = 0; indx_d < velocity_domain_ptr_vec.size(); indx_d++)
    {
        Domain2D *domain_ptr = velocity_domain_ptr_vec[indx_d];
        matrix_fill_velocity_clear(a_mat, b_vec, x_vec, domain_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_domain_ptr_vec.size(); indx_d++)
    {
        Domain2D *domain_ptr = pressure_domain_ptr_vec[indx_d];
        matrix_fill_pressure_clear(a_mat, b_vec, x_vec, domain_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_point_domain_ptr_vec.size(); indx_d++)
    {
        Domain0D *domain_ptr = pressure_point_domain_ptr_vec[indx_d];
        matrix_fill_pressure_point_clear(a_mat, b_vec, x_vec, domain_ptr);
    }
    
    // iterate through each dirichlet boundary
    for (int indx_d = 0; indx_d < velocity_domain_ptr_vec.size(); indx_d++)
    {
         Domain2D *domain_ptr = velocity_domain_ptr_vec[indx_d];
         Scalar2D *velocity_x_ptr = velocity_x_ptr_vec[indx_d];
         Scalar2D *velocity_y_ptr = velocity_y_ptr_vec[indx_d];
         Scalar2D *velocity_z_ptr = velocity_z_ptr_vec[indx_d];
         matrix_fill_velocity(a_mat, b_vec, x_vec, domain_ptr, velocity_x_ptr, velocity_y_ptr, velocity_z_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_domain_ptr_vec.size(); indx_d++)
    {
         Domain2D *domain_ptr = pressure_domain_ptr_vec[indx_d];
         Scalar2D *pressure_ptr = pressure_ptr_vec[indx_d];
         matrix_fill_pressure(a_mat, b_vec, x_vec, domain_ptr, pressure_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_point_domain_ptr_vec.size(); indx_d++)
    {
         Domain0D *domain_ptr = pressure_point_domain_ptr_vec[indx_d];
         Scalar0D *pressure_ptr = pressure_point_ptr_vec[indx_d];
         matrix_fill_pressure_point(a_mat, b_vec, x_vec, domain_ptr, pressure_ptr);
    }

}

void PhysicsSteadyNavierStokes::matrix_fill_domain
(
    std::vector<EigenTriplet> &delta_a_triplet_vec, EigenVector &b_vec, EigenVector &x_vec,
    Domain3D *domain_line_ptr, Domain3D *domain_quad_ptr, IntegralTaylorHood3D *integral_ptr,
    Scalar3D *density_ptr, Scalar3D *viscosity_ptr, Scalar3D *force_x_ptr, Scalar3D *force_y_ptr, Scalar3D *force_z_ptr
)
{

    // offset the starting row depending on the equation
    int offset_nsex = 0;
    int offset_nsey = offset_nsex + velocity_x_ptr->num_point;
    int offset_nsez = offset_nsey + velocity_y_ptr->num_point;
    int offset_cont = offset_nsez + velocity_z_ptr->num_point;

    // fill up matrix with x-component of NSE
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble den_vec = density_ptr->get_neighbor_value(edid);
        VectorDouble visc_vec = viscosity_ptr->get_neighbor_value(edid);
        VectorDouble fcex_vec = force_x_ptr->get_neighbor_value(edid);

        // get velocities of points around element
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble velz_vec = velocity_z_ptr->get_neighbor_value(domain_quad_ptr, edid);

        // get group ID of points
        int egid = domain_quad_ptr->element_edid_to_egid_vec[edid];
        int edid_line = domain_line_ptr->element_egid_to_edid_map[egid];  // local element ID of linear domain
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_quad_ptr, edid);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid_line);

        // matrix indexing
        // matrix row = start_row of test function (physics) + offset + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {

            // velocity x columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate velocity integrals
                double integ_Mi_velx_derv_Mj_x = 0;
                double integ_Mi_vely_derv_Mj_y = 0;
                double integ_Mi_velz_derv_Mj_z = 0;
                for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++)
                {
                    integ_Mi_velx_derv_Mj_x += velx_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_x_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_vely_derv_Mj_y += vely_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_y_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_velz_derv_Mj_z += velz_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_z_quad_vec[edid][indx_i][indx_j][indx_k];
                }
                
                // fill up matrix
                int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
                int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(
                    mat_row, mat_col,
                    den_vec[indx_i] * (integ_Mi_velx_derv_Mj_x + integ_Mi_vely_derv_Mj_y + integ_Mi_velz_derv_Mj_z) +
                    visc_vec[indx_i] * integral_ptr->integral_div_Mi_dot_div_Mj_quad_vec[edid][indx_i][indx_j]
                ));

            }

            // pressure columns
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Mi_derivative_Nj_x_quad_vec[edid][indx_i][indx_j]));
            }
            
        }

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += den_vec[indx_i] * fcex_vec[indx_i] * integral_ptr->integral_Mi_quad_vec[edid][indx_i];
        }

    }

    // fill up matrix with y-component of NSE
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble den_vec = density_ptr->get_neighbor_value(edid);
        VectorDouble visc_vec = viscosity_ptr->get_neighbor_value(edid);
        VectorDouble fcey_vec = force_y_ptr->get_neighbor_value(edid);

        // get velocities of points around element
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble velz_vec = velocity_z_ptr->get_neighbor_value(domain_quad_ptr, edid);

        // get group ID of points
        int egid = domain_quad_ptr->element_edid_to_egid_vec[edid];
        int edid_line = domain_line_ptr->element_egid_to_edid_map[egid];  // local element ID of linear domain
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_quad_ptr, edid);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid_line);

        // matrix indexing
        // matrix row = start_row of test function (physics) + offset + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {

            // velocity y columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate velocity integrals
                double integ_Mi_velx_derv_Mj_x = 0;
                double integ_Mi_vely_derv_Mj_y = 0;
                double integ_Mi_velz_derv_Mj_z = 0;
                for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++)
                {
                    integ_Mi_velx_derv_Mj_x += velx_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_x_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_vely_derv_Mj_y += vely_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_y_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_velz_derv_Mj_z += velz_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_z_quad_vec[edid][indx_i][indx_j][indx_k];
                }
                
                // fill up matrix
                int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
                int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col,
                    den_vec[indx_i] * (integ_Mi_velx_derv_Mj_x + integ_Mi_vely_derv_Mj_y + integ_Mi_velz_derv_Mj_z) +
                    visc_vec[indx_i] * integral_ptr->integral_div_Mi_dot_div_Mj_quad_vec[edid][indx_i][indx_j]
                ));

            }

            // pressure columns
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Mi_derivative_Nj_y_quad_vec[edid][indx_i][indx_j]));
            }
            
        }

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += den_vec[indx_i] * fcey_vec[indx_i] * integral_ptr->integral_Mi_quad_vec[edid][indx_i];
        }

    }

    // fill up matrix with z-component of NSE
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble den_vec = density_ptr->get_neighbor_value(edid);
        VectorDouble visc_vec = viscosity_ptr->get_neighbor_value(edid);
        VectorDouble fcez_vec = force_z_ptr->get_neighbor_value(edid);

        // get velocities of points around element
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble velz_vec = velocity_z_ptr->get_neighbor_value(domain_quad_ptr, edid);

        // get group ID of points
        int egid = domain_quad_ptr->element_edid_to_egid_vec[edid];
        int edid_line = domain_line_ptr->element_egid_to_edid_map[egid];  // local element ID of linear domain
        VectorInt velz_pfid_vec = velocity_z_ptr->get_neighbor_pfid(domain_quad_ptr, edid);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid_line);

        // matrix indexing
        // matrix row = start_row of test function (physics) + offset + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {

            // velocity y columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate velocity integrals
                double integ_Mi_velx_derv_Mj_x = 0;
                double integ_Mi_vely_derv_Mj_y = 0;
                double integ_Mi_velz_derv_Mj_z = 0;
                for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++)
                {
                    integ_Mi_velx_derv_Mj_x += velx_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_x_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_vely_derv_Mj_y += vely_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_y_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_velz_derv_Mj_z += velz_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_z_quad_vec[edid][indx_i][indx_j][indx_k];
                }
                
                // fill up matrix
                int mat_row = start_row + offset_nsez + velz_pfid_vec[indx_i];
                int mat_col = velocity_z_ptr->start_col + velz_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col,
                    den_vec[indx_i] * (integ_Mi_velx_derv_Mj_x + integ_Mi_vely_derv_Mj_y + integ_Mi_velz_derv_Mj_z) +
                    visc_vec[indx_i] * integral_ptr->integral_div_Mi_dot_div_Mj_quad_vec[edid][indx_i][indx_j]
                ));

            }

            // pressure columns
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_nsez + velz_pfid_vec[indx_i];
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Mi_derivative_Nj_z_quad_vec[edid][indx_i][indx_j]));
            }
            
        }

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsez + velz_pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += den_vec[indx_i] * fcez_vec[indx_i] * integral_ptr->integral_Mi_quad_vec[edid][indx_i];
        }

    }

    // fill up matrix with continuity equation
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // get group ID of points
        int egid = domain_line_ptr->element_edid_to_egid_vec[edid];
        int edid_quad = domain_quad_ptr->element_egid_to_edid_map[egid];  // local element ID of quadratic domain
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_quad_ptr, edid_quad);
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_quad_ptr, edid_quad);
        VectorInt velz_pfid_vec = velocity_z_ptr->get_neighbor_pfid(domain_quad_ptr, edid_quad);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid);

        // matrix indexing
        // matrix row = start_row of test function (physics) + offset + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++)
        {

            // velocity x columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
                int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Ni_derivative_Mj_x_line_vec[edid][indx_i][indx_j]));
            }

            // velocity y columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
                int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Ni_derivative_Mj_y_line_vec[edid][indx_i][indx_j]));
            }

            // velocity z columns
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
                int mat_col = velocity_z_ptr->start_col + velz_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Ni_derivative_Mj_z_line_vec[edid][indx_i][indx_j]));
            }

            // pressure columns (for stability)
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];
                delta_a_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, 1e-8 * integral_ptr->integral_Ni_Nj_line_vec[edid][indx_i][indx_j]));
            }

        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_velocity_clear
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_ptr
)
{

    // offset the starting row depending on the equation
    int offset_nsex = 0;
    int offset_nsey = offset_nsex + velocity_x_ptr->num_point;
    int offset_nsez = offset_nsey + velocity_y_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get group ID of points
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_ptr, edid);
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_ptr, edid);
        VectorInt velz_pfid_vec = velocity_z_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows (velocity x)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int velx_mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
            a_mat.row(velx_mat_row) *= 0.;
            b_vec.coeffRef(velx_mat_row) = 0.;
        }

        // clear rows (velocity y)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int vely_mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
            a_mat.row(vely_mat_row) *= 0.;
            b_vec.coeffRef(vely_mat_row) = 0.;
        }

        // clear rows (velocity z)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int velz_mat_row = start_row + offset_nsez + velz_pfid_vec[indx_i];
            a_mat.row(velz_mat_row) *= 0.;
            b_vec.coeffRef(velz_mat_row) = 0.;
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure_clear
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point + velocity_z_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure_point_clear
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain0D *domain_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point + velocity_z_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_velocity
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_ptr,
    Scalar2D *velx_ptr, Scalar2D *vely_ptr, Scalar2D *velz_ptr
)
{

    // offset the starting row depending on the equation
    int offset_nsex = 0;
    int offset_nsey = offset_nsex + velocity_x_ptr->num_point;
    int offset_nsez = offset_nsey + velocity_y_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble velx_vec = velx_ptr->get_neighbor_value(edid);
        VectorDouble vely_vec = vely_ptr->get_neighbor_value(edid);
        VectorDouble velz_vec = velz_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_ptr, edid);
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_ptr, edid);
        VectorInt velz_pfid_vec = velocity_z_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows (velocity x)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
            int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += velx_vec[indx_i];
        }

        // clear rows (velocity y)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
            int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += vely_vec[indx_i];
        }

        // clear rows (velocity z)
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsez + velz_pfid_vec[indx_i];
            int mat_col = velocity_z_ptr->start_col + velz_pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += velz_vec[indx_i];
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_ptr,
    Scalar2D *pres_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point + velocity_z_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble pres_vec = pres_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += pres_vec[indx_i];
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure_point
(
    EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec,
    Domain0D *domain_ptr,
    Scalar0D *pres_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point + velocity_z_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble pres_vec = pres_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += pres_vec[indx_i];
        }

    }

}

}

#endif
