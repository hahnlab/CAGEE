#ifndef PROPORTIONAL_VARIANCE_H
#define PROPORTIONAL_VARIANCE_H

#include <cmath>

/// <summary>
/// Convert values from user space (linear space) to computational space (log space) and back again.
/// See architecture document #12
/// </summary>


namespace proportional_variance
{

    inline double to_user_space(double val)
    {
#ifdef MODEL_GENE_EXPRESSION_LOGS
        return std::exp(val) - LOG_OFFSET;
#else
        return val;
#endif
    }

    inline double to_computational_space(double val)
    {
#ifdef MODEL_GENE_EXPRESSION_LOGS
        return std::log(val + LOG_OFFSET);
#else
        return val;
#endif
    }


}

#endif // PROPORTIONAL_VARIANCE_H
