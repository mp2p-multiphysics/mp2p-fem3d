#ifndef PHYSICSTRANSIENT_DIFFUSION
#define PHYSICSTRANSIENT_DIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_2d.hpp"
#include "domain_3d.hpp"
#include "integral_2d.hpp"
#include "integral_3d.hpp"
#include "physicstransient_base.hpp"
#include "scalar_2d.hpp"
#include "scalar_3d.hpp"
#include "variable_3d.hpp"
#include "variable_group.hpp"

namespace FEM3D
{

class PhysicsTransientDiffusion : public PhysicsTransientBase
{
    /*

    Single-component transient diffusion equation.    
    
    a * du/dt = -div(-b * grad(u)) + c

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_variablegroup : void
        Set variables used in this physics.
    set_domain : void
        Set scalars applied to 2D domains.
    set_boundary_dirichlet : void
        Set a Dirichlet boundary condition along a 1D domain.
    set_boundary_neumann : void
        Set a Neumann boundary condition along a 1D domain.

    */

    public:

    // variables
    VariableGroup* value_ptr;

    // domain objects
    std::vector<Domain3D*> domain_ptr_vec;
    std::vector<Integral3D*> integral_ptr_vec;
    std::vector<Scalar3D*> derivativecoefficient_ptr_vec;
    std::vector<Scalar3D*> diffusioncoefficient_ptr_vec;
    std::vector<Scalar3D*> generationcoefficient_ptr_vec;

    // boundary objects - dirichlet
    std::vector<Domain2D*> dirichlet_domain_ptr_vec;
    std::vector<Scalar2D*> dirichlet_constant_ptr_vec;

    // boundary objects - neumann
    std::vector<Domain2D*> neumann_domain_ptr_vec;
    std::vector<Integral2D*> neumann_integral_ptr_vec;
    std::vector<Scalar2D*> neumann_flux_ptr_vec;

    // vectors of objects to update
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<Scalar3D*> scalar3d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(
        EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt
    );
    void set_variablegroup(VariableGroup &value_in);
    void set_domain(Domain3D &domain_in, Integral3D &integral_in, Scalar3D &derivativecoefficient_in, Scalar3D &diffusioncoefficient_in, Scalar3D &generationcoefficient_in);
    void set_boundary_dirichlet(Domain2D &domain_in, Scalar2D &value_constant_in);
    void set_boundary_neumann(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &value_flux_in);

    // getter and setter functions
    void set_start_row(int start_row_in) {start_row = start_row_in;}
    int get_start_row() {return start_row;}
    std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    std::vector<Scalar3D*> get_scalar3d_ptr_vec() {return scalar3d_ptr_vec;}
    std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsTransientDiffusion() {}

    private:

    void matrix_fill_domain
    (
        std::vector<EigenTriplet> &delta_a_triplet_vec, std::vector<EigenTriplet> &delta_c_triplet_vec, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
        Domain3D *domain_ptr, Integral3D *integral_ptr,
        Scalar3D *derivativecoefficient_ptr, Scalar3D *diffusioncoefficient_ptr, Scalar3D *generationcoefficient_ptr
    );
    void matrix_fill_neumann
    (
        EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
        Domain2D *domain_ptr, Integral2D *integral_ptr,
        Scalar2D *value_flux_ptr
    );
    void matrix_fill_dirichlet_clear
    (
        EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
        Domain2D *domain_ptr
    );
    void matrix_fill_dirichlet
    (
        EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
        Domain2D *domain_ptr,
        Scalar2D *value_constant_ptr
    );

};

void PhysicsTransientDiffusion::set_variablegroup(VariableGroup &value_in)
{
    /*
    
    Set variables used in this physics.

    Arguments
    =========
    value_in : VariableGroup
        u in a * du/dt = -div(-b * grad(u)) + c.
        This will be solved for by the matrix equation.

    Returns
    =======
    (none)
    
    */

    // set variable groups
    value_ptr = &value_in;

    // add to vector of variable groups
    variablegroup_ptr_vec.push_back(&value_in);

}

void PhysicsTransientDiffusion::set_domain(Domain3D &domain_in, Integral3D &integral_in, Scalar3D &derivativecoefficient_in, Scalar3D &diffusioncoefficient_in, Scalar3D &generationcoefficient_in)
{
    /*
    
    Set scalars applied to 2D domains.

    Arguments
    =========
    domain_in : Domain3D
        Domain that this physics applies to.
    integral_in : Integral3D
        Test function integrals over the domains.
    derivativecoefficient_in : Scalar3D
        a in a * du/dt = -div(-b * grad(u)) + c.
    diffusioncoefficient_in : Scalar3D
        b in a * du/dt = -div(-b * grad(u)) + c.
    generationcoefficient_in : Scalar3D
        c in a * du/dt = -div(-b * grad(u)) + c.

    Returns
    =======
    (none)

    */

    // add to vector of domain objects
    domain_ptr_vec.push_back(&domain_in);
    integral_ptr_vec.push_back(&integral_in);
    derivativecoefficient_ptr_vec.push_back(&derivativecoefficient_in);
    diffusioncoefficient_ptr_vec.push_back(&diffusioncoefficient_in);
    generationcoefficient_ptr_vec.push_back(&generationcoefficient_in);

    // add to vector of scalar3d objects
    scalar3d_ptr_vec.push_back(&derivativecoefficient_in);
    scalar3d_ptr_vec.push_back(&diffusioncoefficient_in);
    scalar3d_ptr_vec.push_back(&generationcoefficient_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();
    integral_in.evaluate_integral_Ni_Nj();
    integral_in.evaluate_integral_div_Ni_dot_div_Nj();

}

void PhysicsTransientDiffusion::set_boundary_dirichlet(Domain2D &domain_in, Scalar2D &value_constant_in)
{
    /*
    
    Set a Dirichlet boundary condition along a 1D domain.

    Arguments
    =========
    domain_in : Domain2D
        Domain that this boundary condition applies to.
    value_constant_in : Scalar2D
        Constant value prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    dirichlet_domain_ptr_vec.push_back(&domain_in);
    dirichlet_constant_ptr_vec.push_back(&value_constant_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&value_constant_in);

}

void PhysicsTransientDiffusion::set_boundary_neumann(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &value_flux_in)
{
    /*
    
    Set a Neumann boundary condition along a 1D domain.

    Arguments
    =========
    domain_in : Domain2D
        Domain that this boundary condition applies to.
    integral_in : Integral2D
        Test function integrals over the domains.
    value_flux_in : Scalar2D
        Flux prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    neumann_domain_ptr_vec.push_back(&domain_in);
    neumann_integral_ptr_vec.push_back(&integral_in);
    neumann_flux_ptr_vec.push_back(&value_flux_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&value_flux_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();

}

void PhysicsTransientDiffusion::matrix_fill
(
    EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt
)
{
    /*

    Fill up the matrix equation Ax(t+1) = Cx(t) + d with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : EigenSparseMatrix
        A in Ax(t+1) = Cx(t) + d.
    c_mat : EigenSparseMatrix
        C in Ax(t+1) = Cx(t) + d.
    d_vec : EigenVector
        d in Ax(t+1) = Cx(t) + d.
    x_vec : EigenVector
        x(t+1) in Ax(t+1) = Cx(t) + d.
    x_last_timestep_vec : EigenVector
        x(t) in Ax(t+1) = Cx(t) + d.
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

    // represent matrix as triplets for performance
    std::vector<EigenTriplet> delta_a_triplet_vec;
    std::vector<EigenTriplet> delta_c_triplet_vec;
    int num_equation = a_mat.rows();
    delta_a_triplet_vec.reserve(10*num_equation); // estimated number of entries
    delta_c_triplet_vec.reserve(10*num_equation); // estimated number of entries

    // iterate through each domain
    for (int indx_d = 0; indx_d < domain_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain3D *domain_ptr = domain_ptr_vec[indx_d];
        Integral3D *integral_ptr = integral_ptr_vec[indx_d];
        Scalar3D *derivativecoefficient_ptr = derivativecoefficient_ptr_vec[indx_d];
        Scalar3D *diffusioncoefficient_ptr = diffusioncoefficient_ptr_vec[indx_d];
        Scalar3D *generationcoefficient_ptr = generationcoefficient_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(
            delta_a_triplet_vec, delta_c_triplet_vec, d_vec, x_vec, x_last_timestep_vec, dt,
            domain_ptr, integral_ptr,
            derivativecoefficient_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr
        );

    }

    // convert triplet vector to sparse matrix
    EigenSparseMatrix delta_a_mat(num_equation, num_equation);
    EigenSparseMatrix delta_c_mat(num_equation, num_equation);
    delta_a_mat.setFromTriplets(delta_a_triplet_vec.begin(), delta_a_triplet_vec.end());
    delta_c_mat.setFromTriplets(delta_c_triplet_vec.begin(), delta_c_triplet_vec.end());
    a_mat += delta_a_mat;
    c_mat += delta_c_mat;

    // iterate through each neumann boundary
    for (int indx_d = 0; indx_d < neumann_domain_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain2D *domain_ptr = neumann_domain_ptr_vec[indx_d];
        Integral2D *integral_ptr = neumann_integral_ptr_vec[indx_d];
        Scalar2D *value_flux_ptr = neumann_flux_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_neumann(
            a_mat, c_mat, d_vec, x_vec, x_last_timestep_vec, dt,
            domain_ptr, integral_ptr, value_flux_ptr
        );

    }

    // clear equations with dirichlet boundary conditions
    for (int indx_d = 0; indx_d < dirichlet_domain_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain2D *domain_ptr = dirichlet_domain_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_dirichlet_clear(
            a_mat, c_mat, d_vec, x_vec, x_last_timestep_vec, dt,
            domain_ptr
        );

    }

    // iterate through each dirichlet boundary
    for (int indx_d = 0; indx_d < dirichlet_domain_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain2D *domain_ptr = dirichlet_domain_ptr_vec[indx_d];
        Scalar2D *value_constant_ptr = dirichlet_constant_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_dirichlet(
            a_mat, c_mat, d_vec, x_vec, x_last_timestep_vec, dt,
            domain_ptr, value_constant_ptr
        );

    }

}

void PhysicsTransientDiffusion::matrix_fill_domain
(
    std::vector<EigenTriplet> &delta_a_triplet_vec, std::vector<EigenTriplet> &delta_c_triplet_vec, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
    Domain3D *domain_ptr, Integral3D *integral_ptr,
    Scalar3D *derivativecoefficient_ptr, Scalar3D *diffusioncoefficient_ptr, Scalar3D *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble dervcoeff_vec = derivativecoefficient_ptr->get_neighbor_value(edid);
        VectorDouble diffcoeff_vec = diffusioncoefficient_ptr->get_neighbor_value(edid);
        VectorDouble gencoeff_vec = generationcoefficient_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // matrix indexing
        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat and c_mat coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){
        for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_ptr->start_col + pfid_vec[indx_j];
            delta_a_triplet_vec.push_back(EigenTriplet(
                mat_row, mat_col,
                (dervcoeff_vec[indx_i]/dt) * integral_ptr->integral_Ni_Nj_vec[edid][indx_i][indx_j] +
                diffcoeff_vec[indx_i] * integral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j]
            ));
            delta_c_triplet_vec.push_back(EigenTriplet(mat_row, mat_col, (dervcoeff_vec[indx_i]/dt) * integral_ptr->integral_Ni_Nj_vec[edid][indx_i][indx_j]));
        }}

        // calculate d_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            d_vec.coeffRef(mat_row) += gencoeff_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_neumann
(
    EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
    Domain2D *domain_ptr, Integral2D *integral_ptr,
    Scalar2D *value_flux_ptr
)
{

    // iterate for each flux boundary element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {
        
        // get scalar values of points around element
        VectorDouble value_flux_vec = value_flux_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            d_vec.coeffRef(mat_row) += value_flux_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_dirichlet_clear
(
    EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
    Domain2D *domain_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            a_mat.row(mat_row) *= 0.;
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;
        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_dirichlet
(
    EigenSparseMatrix &a_mat, EigenSparseMatrix &c_mat, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt,
    Domain2D *domain_ptr,
    Scalar2D *value_constant_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble value_constant_vec = value_constant_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_ptr->start_col + pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            d_vec.coeffRef(mat_row) += value_constant_vec[indx_i];
        }

    }
    
}

}

#endif
