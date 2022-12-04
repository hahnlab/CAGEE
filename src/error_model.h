#ifndef ERROR_MODEL_H
#define ERROR_MODEL_H

#include <vector>
#include <string>
#include <map>

class error_model {
private:
    size_t _max_family_size;
    std::vector<int> _deviations;
    
     //!< Deviations from the true gene family (e.g., -1 0 1)
  
    std::vector<std::vector<double> > _error_dists; //!< Each vector element will be a gene family size; the vector of doubles inside (e.g., 0.1 0.8 0.1) will be the probs of deviating from the true value
    
public:
    error_model();
    //! Set max family size for which deviations apply
    void set_max_family_size(size_t max_cnt);

    //! Calculate normal Distribution
    double Normal(int x, int mu, double sigma);

    //! Set deviations
    void set_deviations(std::map<int,double> deviations);

    //! Sum of vectors
    double Sum(std::map<int,double> dic);


    //!generate matrix
    std::map<int,double> generate_matrix(double a,double r,double mu,int upper_bound,double ux);

    //! ReNormalize the vectors
    std::map<int,double> reNormalize(std::map<int,double>  dic);

    //! Set deviation probability vector for a certain family size
    void set_probabilities(double a, double r, size_t mu, double upperbound,double ux);

    //! Get deviation probability vector for a certain family size
     std::vector<double> get_probs(size_t mu) const;

    size_t n_deviations() const {
        return _deviations.size();
    }

    size_t get_max_family_size() const {
        return _error_dists.size();
    }
    std::vector<double> get_epsilons() const;
    //void replace_epsilons(std::map<double, double>* new_epsilons);
    void update_single_epsilon(double new_epsilon);

    friend void write_error_model_file(std::ostream& ost, error_model& errormodel);
};

#endif

