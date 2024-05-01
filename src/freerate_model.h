#ifndef FREERATE_MODEL_H
#define FREERATE_MODEL_H

#include "core.h"

class freerate_model : public model {

public:
    //! Computation or estimation constructor
    freerate_model();
    
    virtual ~freerate_model() {}

    virtual double infer_transcript_likelihoods(const user_data& ud, const sigma_squared*p_sigma) override;

    virtual sigma_optimizer_scorer* get_sigma_optimizer(const user_data& data, const std::vector<std::string>& sample_groups) override;

    virtual reconstruction* reconstruct_ancestral_states(const user_data& ud, matrix_cache *p_calc);
};

#endif
