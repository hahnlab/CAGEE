#include <algorithm>
#include <iomanip>
#include <sstream>

#include "lambda.h"
#include "matrix_cache.h"
#include "clade.h"

using namespace std;

void lambda::update(const double* values)
{
    std::copy(values, values + _lambdas.size(), _lambdas.begin());
}

std::string lambda::to_string() const
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

bool lambda::is_valid() const
{
    return std::none_of(_lambdas.begin(), _lambdas.end(), [](double d) { return d < 0; });
}

double lambda::get_value_for_clade(const clade *c) const {
    if (count() == 1)
        return _lambdas[0];

    int index = _node_name_to_lambda_index.at(c->get_taxon_name());
    return _lambdas[index];
}

/* END: Holding lambda values and specifying how likelihood is computed depending on the number of different lambdas */

