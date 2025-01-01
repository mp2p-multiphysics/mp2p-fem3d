#ifndef PHYSICSSTEADY_DIFFUSION
#define PHYSICSSTEADY_DIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_2d.hpp"
#include "domain_3d.hpp"
#include "integral_2d.hpp"
#include "integral_3d.hpp"
#include "physicssteady_base.hpp"
#include "scalar_2d.hpp"
#include "scalar_3d.hpp"
#include "variable_3d.hpp"
#include "variable_group.hpp"

namespace FEM3D
{

class PhysicsSteadyDiffusion : public PhysicsSteadyBase
{
    /*

    Single-component steady-state diffusion equation.    
    
    0 = -div(-b * grad(u)) + c

    Variables
    =========
    value_in : VariableGroup
        u in 0 = -div(-b * grad(u)) + c.
        This will be solved for by the matrix equation.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_variablegroup : void
        Set variables used in this physics.
    set_domain : void
        Set scalars applied to 3D domains.
    set_boundary_dirichlet : void
        Set a Dirichlet boundary condition along a 2D domain.
    set_boundary_neumann : void
        Set a Neumann boundary condition along a 2D domain.

    */

    public:

    // variables
    VariableGroup* value_ptr;

    // vectors indicating dirichlet BCs
    // [pfid] -> true if dirichlet bc
    std::vector<bool> is_value_dirichlet_vec;

    // domain objects
    std::vector<Domain3D*> domain_ptr_vec;
    std::vector<Integral3D*> integral_ptr_vec;
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
    void matrix_fill(EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec);
    void set_domain(Domain3D &domain_in, Integral3D &integral_in, Scalar3D &diffusioncoefficient_in, Scalar3D &generationcoefficient_in);
    void set_boundary_dirichlet(Domain2D &domain_in, Scalar2D &value_constant_in);
    void set_boundary_neumann(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &value_flux_in);

    // getter and setter functions
    void set_start_row(int start_row_in) {start_row = start_row_in;}
    int get_start_row() {return start_row;}
    std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    std::vector<Scalar3D*> get_scalar3d_ptr_vec() {return scalar3d_ptr_vec;}
    std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsSteadyDiffusion() {}

    // constructor
    PhysicsSteadyDiffusion(VariableGroup &value_in)
    {

        // set variable groups
        value_ptr = &value_in;

        // add to vector of variable groups
        variablegroup_ptr_vec.push_back(&value_in);

        // vector indicating if dirichlet bc is applied
        is_value_dirichlet_vec = std::vector<bool> (value_in.num_point, false);

    }

    private:

    void matrix_fill_domain
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain3D *domain_ptr, Integral3D *integral_ptr,
        Scalar3D *diffusioncoefficient_ptr, Scalar3D *generationcoefficient_ptr
    );
    void matrix_fill_neumann
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr, Integral2D *integral_ptr,
        Scalar2D *value_flux_ptr
    );
    void matrix_fill_dirichlet
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_ptr,
        Scalar2D *value_constant_ptr
    );

};

void PhysicsSteadyDiffusion::set_domain(Domain3D &domain_in, Integral3D &integral_in, Scalar3D &diffusioncoefficient_in, Scalar3D &generationcoefficient_in)
{
    /*
    
    Set scalars applied to 3D domains.

    Arguments
    =========
    domain_in : Domain3D
        Domain that this physics applies to.
    integral_in : Integral3D
        Test function integrals over the domains.
    diffusioncoefficient_in : Scalar3D
        b in 0 = -div(-b * grad(u)) + c.
    generationcoefficient_in : Scalar3D
        c in 0 = -div(-b * grad(u)) + c.

    Returns
    =======
    (none)

    */

    // add to vector of domain objects
    domain_ptr_vec.push_back(&domain_in);
    integral_ptr_vec.push_back(&integral_in);
    diffusioncoefficient_ptr_vec.push_back(&diffusioncoefficient_in);
    generationcoefficient_ptr_vec.push_back(&generationcoefficient_in);

    // add to vector of scalar2d objects
    scalar3d_ptr_vec.push_back(&diffusioncoefficient_in);
    scalar3d_ptr_vec.push_back(&generationcoefficient_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();
    integral_in.evaluate_integral_div_Ni_dot_div_Nj();

}

void PhysicsSteadyDiffusion::set_boundary_dirichlet(Domain2D &domain_in, Scalar2D &value_constant_in)
{
    /*
    
    Set a Dirichlet boundary condition along a 2D domain.

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

    // add to vector of scalar1d objects
    scalar2d_ptr_vec.push_back(&value_constant_in);

    // specify points with dirichlet BC
    for (auto pgid : domain_in.point_pdid_to_pgid_vec)
    {
        int pfid = value_ptr->point_pgid_to_pfid_map[pgid];
        is_value_dirichlet_vec[pfid] = true;
    }

}

void PhysicsSteadyDiffusion::set_boundary_neumann(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &value_flux_in)
{
    /*
    
    Set a Neumann boundary condition along a 2D domain.

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

    // add to vector of scalar1d objects
    scalar2d_ptr_vec.push_back(&value_flux_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();

}

void PhysicsSteadyDiffusion::matrix_fill
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec
)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_trivec : EigenTripletVector
        A in Ax = b.
    b_vec : EigenVector
        b in Ax = b.
    x_vec : EigenVector
        x in Ax = b.

    Returns
    =======
    (none)

    */

    // iterate through each domain
    for (int indx_d = 0; indx_d < domain_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain3D *domain_ptr = domain_ptr_vec[indx_d];
        Integral3D *integral_ptr = integral_ptr_vec[indx_d];
        Scalar3D *diffusioncoefficient_ptr = diffusioncoefficient_ptr_vec[indx_d];
        Scalar3D *generationcoefficient_ptr = generationcoefficient_ptr_vec[indx_d];
        
        // fill up matrix with domain equations
        matrix_fill_domain(a_trivec, b_vec, x_vec, domain_ptr, integral_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr);

    }

    // iterate through each neumann boundary
    for (int indx_d = 0; indx_d < neumann_domain_ptr_vec.size(); indx_d++)
    {
 
        // subset domain objects
        Domain2D *domain_ptr = neumann_domain_ptr_vec[indx_d];
        Integral2D *integral_ptr = neumann_integral_ptr_vec[indx_d];
        Scalar2D *value_flux_ptr = neumann_flux_ptr_vec[indx_d];
 
        // fill up matrix with boundary conditions
        matrix_fill_neumann(a_trivec, b_vec, x_vec, domain_ptr, integral_ptr, value_flux_ptr);
 
    }

    // iterate through each dirichlet boundary
    for (int indx_d = 0; indx_d < dirichlet_domain_ptr_vec.size(); indx_d++)
    {
 
        // subset domain objects
        Domain2D *domain_ptr = dirichlet_domain_ptr_vec[indx_d];
        Scalar2D *value_constant_ptr = dirichlet_constant_ptr_vec[indx_d];
 
        // fill up matrix with boundary conditions
        matrix_fill_dirichlet(a_trivec, b_vec, x_vec, domain_ptr, value_constant_ptr);
 
    }

}

void PhysicsSteadyDiffusion::matrix_fill_domain
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain3D *domain_ptr, Integral3D *integral_ptr,
    Scalar3D *diffusioncoefficient_ptr, Scalar3D *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble diffcoeff_vec = diffusioncoefficient_ptr->get_neighbor_value(edid);
        VectorDouble gencoeff_vec = generationcoefficient_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // matrix indexing
        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // iterate through test functions
        // associated with matrix row
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {

            // skip if matrix row has dirichlet condition
            if (is_value_dirichlet_vec[pfid_vec[indx_i]])
            {
                continue;
            }
                
            // calculate matrix row
            int mat_row = start_row + pfid_vec[indx_i];

            // iterate through trial functions
            // associated with matrix column
            for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = value_ptr->start_col + pfid_vec[indx_j];

                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col, diffcoeff_vec[indx_i] * integral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j]));

            }

            // append to b_vec
            b_vec.coeffRef(mat_row) += gencoeff_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];

        }

    }

}

void PhysicsSteadyDiffusion::matrix_fill_neumann
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_ptr, Integral2D *integral_ptr,
    Scalar2D *value_flux_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble value_flux_vec = value_flux_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // iterate through test functions
        // associated with matrix row
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {

            // skip if matrix row has dirichlet condition
            if (is_value_dirichlet_vec[pfid_vec[indx_i]])
            {
                continue;
            }
                
            // calculate matrix row
            int mat_row = start_row + pfid_vec[indx_i];

            // append to b_vec
            b_vec.coeffRef(mat_row) += value_flux_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];

        }

    }

}

void PhysicsSteadyDiffusion::matrix_fill_dirichlet
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
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

        // apply dirichlet BC
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_ptr->start_col + pfid_vec[indx_i];
            a_trivec.push_back(EigenTriplet(mat_row, mat_col, 1.));
            b_vec.coeffRef(mat_row) += value_constant_vec[indx_i];
        }

    }

}

}

#endif
