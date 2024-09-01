#ifndef SCALAR_HEX8
#define SCALAR_HEX8
#include <vector>
#include "Eigen/Eigen"
#include "grid_hex8.hpp"

class ScalarHex8Class
{

    public:
    
    // variables
    GridHex8Struct gh8s;
    std::vector<double> scalar_vec;

    // functions
    void set_from_constant(double value);
    void set_from_vector(std::vector<double> value_vec, int start_id);
    void set_from_evector(Eigen::VectorXd value_evec, int start_id);

    // constructor
    ScalarHex8Class(GridHex8Struct &gh8s_in)
    {
        gh8s = gh8s_in;
    }

    ScalarHex8Class(GridHex8Struct &gh8s_in, double value)
    {
        gh8s = gh8s_in;
        set_from_constant(value);
    }

    ScalarHex8Class(GridHex8Struct &gh8s_in, std::vector<double> value_vec, int start_id)
    {
        gh8s = gh8s_in;
        set_from_vector(value_vec, start_id);
    }

    ScalarHex8Class(GridHex8Struct &gh8s_in, Eigen::VectorXd value_evec, int start_id)
    {
        gh8s = gh8s_in;
        set_from_evector(value_evec, start_id);
    }

};

void ScalarHex8Class::set_from_constant(double value)
{

    // set each scalar_vec value to the constant
    for (int n = 0; n < gh8s.num_point; n++)
    {
        scalar_vec.push_back(value);
    }

}

void ScalarHex8Class::set_from_vector(std::vector<double> value_vec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gh8s.num_point; n++)
    {
        scalar_vec.push_back(value_vec[start_id + n]);
    }

}

void ScalarHex8Class::set_from_evector(Eigen::VectorXd value_evec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gh8s.num_point; n++)
    {
        scalar_vec.push_back(value_evec.coeffRef(start_id + n));
    }

}

#endif