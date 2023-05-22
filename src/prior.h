#ifndef PRIOR_H
#define PRIOR_H

#include <string>

class prior
{
    std::string distribution;
    double param1 = 0.0, param2 = 0.0;
public:
    prior(std::string dist, double p1, double p2) : distribution(dist), param1(p1), param2(p2)
    {

    }
    prior() {}
    double pdf(double value) const;
};

#endif
