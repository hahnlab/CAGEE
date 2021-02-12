#include <vector>
#include <iostream>
#include <valarray>
#include <fstream>
#include <assert.h>
#include <numeric>

#include "easylogging++.h"

#include "core.h"
#include "user_data.h"
#include "matrix_cache.h"
#include "bounded_brownian_motion_model.h"
#include "error_model.h"

std::vector<model *> build_models(const input_parameters& user_input, user_data& user_data) {

    std::vector<gene_family>* p_gene_families = &user_data.gene_families;

    model *p_model = new bounded_brownian_motion_model(user_data.p_lambda, user_data.p_tree, p_gene_families, user_data.max_family_size, user_data.max_root_family_size, nullptr);

    return std::vector<model *>{p_model};
}

std::ostream& operator<<(std::ostream& ost, const family_info_stash& r)
{
    ost << r.family_id << "\t" << r.lambda_multiplier << "\t" << r.category_likelihood << "\t" << r.family_likelihood;
    ost << "\t" << r.posterior_probability << "\t" << (r.significant ? "*" : "N/S");
    return ost;
}

model::model(lambda* p_lambda,
    const clade *p_tree,
    const vector<gene_family> *p_gene_families,
    int max_family_size,
    int max_root_family_size,
    error_model *p_error_model) :
    _ost(cout), _p_lambda(p_lambda), _p_tree(p_tree), _p_gene_families(p_gene_families), _max_family_size(max_family_size),
    _max_root_family_size(max_root_family_size), _p_error_model(p_error_model) 
{
    if (_p_gene_families)
        references = build_reference_list(*_p_gene_families);
}

std::size_t model::get_gene_family_count() const {
    return _p_gene_families->size();
}

void model::initialize_lambda(clade *p_lambda_tree)
{
    lambda *p_lambda = NULL;
    if (p_lambda_tree != NULL)
    {
        std::set<int> unique_lambdas;
        auto fn = [&unique_lambdas](const clade *p_node) { unique_lambdas.insert(p_node->get_lambda_index()); };
        p_lambda_tree->apply_prefix_order(fn);
        auto node_name_to_lambda_index = p_lambda_tree->get_lambda_index_map();
        p_lambda = new multiple_lambda(node_name_to_lambda_index, std::vector<double>(unique_lambdas.size()));
        LOG(INFO) << "Searching for " << unique_lambdas.size() << " lambdas" << endl;
    }
    else
    {
        p_lambda = new single_lambda(0.0);
    }

    _p_lambda = p_lambda;
}

lambda* model::get_simulation_lambda()
{
    return _p_lambda->clone();
}

void model::write_error_model(std::ostream& ost) const
{
    auto em = _p_error_model;
    if (!em)
    {
        em = new error_model();
        em->set_probabilities(_max_family_size, { 0, 1, 0 });
    }
    write_error_model_file(ost, *em);
}

//! Computes likelihoods for the given tree and a single family. Uses a lambda value based on the provided lambda
/// and a given multiplier. Works by calling \ref compute_node_probability on all nodes of the tree
/// using the species counts for the family. 
/// \returns a vector of probabilities for gene counts at the root of the tree 
std::vector<double> inference_prune(const gene_family& gf, matrix_cache& calc, const lambda *p_lambda, const error_model* p_error_model, const clade *p_tree, double lambda_multiplier, int max_root_family_size, int max_family_size)
{
    unique_ptr<lambda> multiplier(p_lambda->multiply(lambda_multiplier));
    clademap<std::vector<double>> probabilities;
    auto init_func = [&](const clade* node) { probabilities[node].resize(node->is_root() ? max_root_family_size : max_family_size + 1); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), init_func);

    auto compute_func = [&](const clade *c) { compute_node_probability(c, gf, p_error_model, probabilities, pair<int, int>(1, max_root_family_size), max_family_size, multiplier.get(), calc); };
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), compute_func);

    return probabilities.at(p_tree); // likelihood of the whole tree = multiplication of likelihood of all nodes
}

bool branch_probabilities::contains(const gene_family& fam) const { 
    return _probabilities.find(fam.id()) != _probabilities.end(); 
}

branch_probabilities::branch_probability branch_probabilities::at(const gene_family& fam, const clade* c) const {
    return _probabilities.at(fam.id()).at(c);
}

void branch_probabilities::set(const gene_family& fam, const clade* c, branch_probability p)
{
    _probabilities[fam.id()][c] = p;
}
