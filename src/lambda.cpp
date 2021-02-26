#include <algorithm>
#include <iomanip>
#include <sstream>

#include "lambda.h"
#include "matrix_cache.h"
#include "clade.h"

using namespace std;

vector<double> get_lambda_values(const lambda* p_lambda)
{
    vector<double> lambdas;
    auto sl = dynamic_cast<const single_lambda*>(p_lambda);
    if (sl)
    {
        lambdas.push_back(sl->get_single_lambda());
    }
    else
    {
        auto ml = dynamic_cast<const multiple_lambda*>(p_lambda);
        lambdas = ml->get_lambdas();
    }
    return lambdas;
}


std::string single_lambda::to_string() const
{
    ostringstream ost;
    ost << setw(15) << setprecision(14) << _lambda;
    return ost.str();
}

void multiple_lambda::update(const double* values)
{
    std::copy(values, values + _lambdas.size(), _lambdas.begin());
}

std::string multiple_lambda::to_string() const
{
    ostringstream ost;
    ost << setw(15) << setprecision(14);
    for (size_t i = 0; i < _lambdas.size(); ++i)
    {
        ost << _lambdas[i];
        if (i != _lambdas.size() - 1) ost << ", ";
    }
    return ost.str();
}

bool multiple_lambda::is_valid()
{
    return std::none_of(_lambdas.begin(), _lambdas.end(), [](double d) { return d < 0; });
}

double multiple_lambda::get_value_for_clade(const clade *c) const {
    int index = _node_name_to_lambda_index.at(c->get_taxon_name());
    return _lambdas[index];
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

