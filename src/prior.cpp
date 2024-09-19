#include <random>

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/normal.hpp>

#include "doctest.h"
#include "easylogging++.h"

#include "prior.h"
#include "matrix_cache.h"
#include "clade.h"
#include "gene_transcript.h"

using namespace std;

prior::prior(std::string dist, double p1, double p2) : distribution(dist), param1(p1), param2(p2)
{
    auto supported_distributions = {"gamma", "fisher", "uniform", "normal"};
    if (std::find(supported_distributions.begin(), supported_distributions.end(), dist) == supported_distributions.end())
        throw std::domain_error("Prior must be given in the form dist:param1:param2");
}


double prior::pdf(double value) const
{
    if (distribution == "gamma")
    {
        boost::math::gamma_distribution<double> d(param1, param2);
        return boost::math::pdf(d, value);
    }
    else if (distribution == "fisher")
    {
        boost::math::fisher_f_distribution<double> f(param1, param2);
        return boost::math::pdf(f, value);
    }
    else if (distribution == "uniform")
    {
        boost::math::uniform_distribution<double> f(param1, param2);
        return boost::math::pdf(f, value);
    }
    else if (distribution == "normal")
    {
        boost::math::normal_distribution<double> f(param1, param2);
        return boost::math::pdf(f, value);
    }
    else
    {
        throw std::domain_error("Unknown probability distribution type");
    }
}

vector<double> get_priors(const matrix_cache& calc, boundaries bounds, const prior *p_prior)
{
    vector<double> priors(calc.create_vector().size());
    for (size_t j = 0; j < priors.size(); ++j) {
        double x = (double(j) + 0.5) * double(bounds.second) / (priors.size() - 1);

        priors[j] = computational_space_prior(x, p_prior);
    }

    return priors;
}

double computational_space_prior(double val, const prior *p_prior)
{
#ifdef MODEL_GENE_EXPRESSION_LOGS
    return exp(val) * p_prior->pdf(exp(val));
#else
    return p_prior->pdf(val);
#endif

}

double compute_prior_likelihood(const vector<double>& partial_likelihood, const vector<double>& priors)
{
    std::vector<double> full(partial_likelihood.size());
    std::transform(partial_likelihood.begin(), partial_likelihood.end(), priors.begin(), full.begin(), std::multiplies<double>());
    std::transform(full.begin(), full.end(), full.begin(), [](double d) {
        return isnan(d) ? -numeric_limits<double>::infinity() : d;
        });

#ifdef USE_MAX_PROBABILITY
    double likelihood = *max_element(full.begin(), full.end()); // get max (CAFE's approach)
#else
    double likelihood = accumulate(full.begin(), full.end(), 0.0, [](double a, double b) { return isinf(b) ? a : a+b; }); // sum over all sizes (Felsenstein's approach)
#endif
    return likelihood;
}

prior estimate_distribution(string dist, const std::vector<double>& data) {
    if (data.size() < 2)    // if there is only one data point, return a weak prior
    {
        if (dist == "gamma")
            return prior("gamma", 0.375, 1600);
        else if (dist == "normal")
            return prior("normal", 1.0, 1.0);
    }

    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();

    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    double variance = sq_sum / data.size() - mean * mean;

    if (dist == "gamma")
    {
        double shape = mean * mean / variance; // k
        double scale = variance / mean; // theta
        return prior("gamma", shape, scale);
    }
    else if (dist == "normal")
    {
        return prior("normal", mean, sqrt(variance));
    }
    else
    {
        throw std::domain_error("Unknown distribution type");
    }   
}

/// @brief  Estimate the distribution of a vector of values
/// @param p_tree 
/// @param gene_transcripts 
/// @param ratios 
/// @return 
clademap<prior> compute_tree_priors(const clade* p_tree, const vector<gene_transcript>& gene_transcripts, bool ratios)
{
    clademap<vector<double>> averages;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&](const clade* c) {
        if (!c->is_leaf())
        {
            averages[c] = vector<double>(gene_transcripts.size());
            int child_count = 0;
            c->apply_to_descendants([&](const clade* d) 
            { 
                if (averages.find(d) != averages.end())
                {
                    // child has already been processed, so add its (already averaged) values to the current node
                    transform(averages[d].begin(), averages[d].end(), averages[c].begin(), averages[c].begin(), std::plus<double>());
                }
                else
                {
                    // child is a leaf, so add its transcript values to the current node
                    transform(gene_transcripts.begin(), gene_transcripts.end(), averages[c].begin(), averages[c].begin(), [d, ratios](const gene_transcript &gf, double val) {
                        double newval = gf.get_expression_value(d->get_taxon_name());
                        return val + (ratios ? newval : exp(newval));
                    });
                }
                child_count++;  
            });
            transform(averages[c].begin(), averages[c].end(), averages[c].begin(), [child_count](double val) { return val / child_count; });  
        }
    });
    clademap<prior> result;
    for(auto& a : averages)
    {
        string dist = ratios ? "normal" : "gamma";
        result[a.first] = estimate_distribution(dist, a.second);
    }
    return result;
}


std::ostream& operator<<(std::ostream& ost, const prior& p)
{
    ost << p.distribution << ":" << p.param1 << ":" << p.param2;
    return ost;
}

TEST_CASE("prior returns correct pdf values for gamma")
{
    prior p("gamma", 0.375, 1600);
    CHECK_EQ(doctest::Approx(0.0133233), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.00719489), p.pdf(8));
    CHECK_EQ(doctest::Approx(0.00505239), p.pdf(14));
}

TEST_CASE("prior returns correct pdf values for fisher")
{
    prior p("fisher", 0.75, 0.75);
    CHECK_EQ(doctest::Approx(0.0388044), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.0114423), p.pdf(8));
    CHECK_EQ(doctest::Approx(0.0054983), p.pdf(14));
}

TEST_CASE("prior returns correct pdf values for uniform")
{
    prior p("uniform", 0, 10);
    CHECK_EQ(doctest::Approx(0.1), p.pdf(3));
    CHECK_EQ(doctest::Approx(0.1), p.pdf(8));
    CHECK_EQ(doctest::Approx(0), p.pdf(14));
}

TEST_CASE("prior throws on unknown")
{
    CHECK_THROWS_AS(prior("", 0,0), std::domain_error);
    prior p;
    CHECK_THROWS_AS(p.pdf(3), std::domain_error);
}

TEST_CASE("get_priors")
{
    matrix_cache calc;
    auto p_prior = new prior("gamma", 1.0, 1600);
    auto bounds = boundaries(0, 20);
    auto priors = get_priors(calc, bounds, p_prior);
    REQUIRE_EQ(200, priors.size());
    CHECK_EQ(doctest::Approx(-7.32816), log(priors[0]));
    CHECK_EQ(doctest::Approx(-1.02372), log(priors[75]));
    CHECK_EQ(doctest::Approx(-182.551), log(priors[125]));
    CHECK_EQ(0, priors[199]);
}

TEST_CASE("compute_prior_likelihood combines prior and inference correctly")
{
    vector<double> inf{ 0.1, 0.2, 0.3};

    vector<double> priors({ 1.43078e-15,    2.5363e-23,  5.65526e-35 });
    double actual = log(compute_prior_likelihood(inf, priors));

#ifdef USE_MAX_PROBABILITY
    CHECK_EQ(doctest::Approx(-35.7683), actual);
#else
    CHECK_EQ(doctest::Approx(-36.4831), actual);
#endif
}

TEST_CASE("prior to string")
{
    prior p("gamma", 0.375, 1600);
    ostringstream ost;
    ost << p;
    CHECK_EQ("gamma:0.375:1600", ost.str());
}

TEST_CASE("estimate_distribution returns correct gamma prior")
{
    vector<double> data{ 0.1, 0.2, 0.3, 0.4, 0.5 };
    auto p = estimate_distribution("gamma", data);
    ostringstream ost;
    ost << p;
    CHECK_EQ("gamma:4.5:0.0666667", ost.str());
}

TEST_CASE("estimate_distribution returns correct normal prior")
{
    vector<double> data{ 0.1, 0.2, 0.3, 0.4, 0.5 };
    auto p = estimate_distribution("normal", data);
    ostringstream ost;
    ost << p;
    CHECK_EQ("normal:0.3:0.141421", ost.str());
}

#define CHECK_STREAM_CONTAINS(x,y) CHECK_MESSAGE(x.str().find(y) != std::string::npos, x.str())

TEST_CASE("node_priors")
{
    unique_ptr<clade> t(parse_newick("((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;"));

    vector<gene_transcript> gene_transcripts;
    gene_transcripts.push_back(gene_transcript("TestFamily1", "", ""));
    gene_transcripts.push_back(gene_transcript("TestFamily2", "", ""));
    gene_transcripts[0].set_expression_value("A", 0.5);
    gene_transcripts[0].set_expression_value("B", 1);
    gene_transcripts[0].set_expression_value("C", 1.5);
    gene_transcripts[0].set_expression_value("D", 2);
    gene_transcripts[0].set_expression_value("E", 3);
    gene_transcripts[1].set_expression_value("A", 8);
    gene_transcripts[1].set_expression_value("B", 4);
    gene_transcripts[1].set_expression_value("C", 4.5);
    gene_transcripts[1].set_expression_value("D", 3.5);
    gene_transcripts[1].set_expression_value("E", 2.5);

    auto result = compute_tree_priors(t.get(), gene_transcripts, false);
    CHECK_EQ(4, result.size());
    ostringstream ost;
    for(auto& a : result)
    {
        ost << a.first->get_ape_index() << ":" << a.second << "|";
    }   
    CHECK_STREAM_CONTAINS(ost, "6:gamma:1.08613:194.18");
    CHECK_STREAM_CONTAINS(ost, "7:gamma:16.6708:1.09132");
    CHECK_STREAM_CONTAINS(ost, "8:gamma:1.01672:396.977");
    CHECK_STREAM_CONTAINS(ost, "9:gamma:1.00577:755.62");
}
