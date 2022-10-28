#ifndef REPLICATE_MODEL_H
#define REPLICATE_MODEL_H

#include <map>
#include <string>

class replicate_model
{
public:
    // replicate id to species id
    std::map<std::string, std::string> _replicates;

    void apply() const;
};


















#endif