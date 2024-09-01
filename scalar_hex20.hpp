#ifndef SCALAR_HEX20
#define SCALAR_HEX20
#include <vector>
#include "Eigen/Eigen"
#include "grid_hex20.hpp"

class ScalarHex20Class
{

    public:
    
    // variables
    GridHex20Struct gh20s;
    std::vector<double> scalar_vec;

    // functions
    void set_from_constant(double value);
    void set_from_vector(std::vector<double> value_vec, int start_id);
    void set_from_evector(Eigen::VectorXd value_evec, int start_id);

    // constructor
    ScalarHex20Class(GridHex20Struct &gh20s_in)
    {
        gh20s = gh20s_in;
    }

    ScalarHex20Class(GridHex20Struct &gh20s_in, double value)
    {
        gh20s = gh20s_in;
        set_from_constant(value);
    }

    ScalarHex20Class(GridHex20Struct &gh20s_in, std::vector<double> value_vec, int start_id)
    {
        gh20s = gh20s_in;
        set_from_vector(value_vec, start_id);
    }

    ScalarHex20Class(GridHex20Struct &gh20s_in, Eigen::VectorXd value_evec, int start_id)
    {
        gh20s = gh20s_in;
        set_from_evector(value_evec, start_id);
    }

};

void ScalarHex20Class::set_from_constant(double value)
{

    // set each scalar_vec value to the constant
    for (int n = 0; n < gh20s.num_point; n++)
    {
        scalar_vec.push_back(value);
    }

}

void ScalarHex20Class::set_from_vector(std::vector<double> value_vec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gh20s.num_point; n++)
    {
        scalar_vec.push_back(value_vec[start_id + n]);
    }

}

void ScalarHex20Class::set_from_evector(Eigen::VectorXd value_evec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gh20s.num_point; n++)
    {
        scalar_vec.push_back(value_evec.coeffRef(start_id + n));
    }

}

#endif