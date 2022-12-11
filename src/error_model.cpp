#include <numeric>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "error_model.h"
#include <iostream>

using namespace std;

error_model::error_model()
{
    _deviations = {};
}

void error_model::set_max_family_size(size_t max_cnt) {
    _max_family_size = max_cnt;
}

void error_model::set_deviations(std::map<double,double> deviations) {
    std::vector<int> deviations1;
    for(auto i = deviations.begin(); i!=deviations.end();++i)
    {   _deviations.push_back(i->first);
    }


}

inline bool is_nearly_equal(double x, double y)
{
    const double epsilon = 0.01;
    return std::abs(x - y) <= epsilon * std::abs(x);
}

double error_model::Normal(double x, double mu, double sigma)
{
double deno,numo,density;
deno = pow((2*M_PI*sigma*sigma),-0.5);
numo = exp(-0.5*pow(((x-mu)/sigma),2));
density = numo*deno;
return density;

}

double error_model::Sum(std::map<double,double> dic)
{
    double sum_=0;
     for(auto i = dic.begin(); i!=dic.end();++i)
{    sum_=sum_+i->second;}

return sum_;
}


std::map<double,double> error_model::reNormalize(std::map<double,double>  dic){
	double sum_values = Sum(dic);
    for(auto i = dic.begin(); i!=dic.end();++i)
{   i->second= i->second/sum_values;
}

	return dic;
}



std::map<double,double> error_model::generate_matrix(double a,double r,double mu,double upper_bound){
	double sigma= exp(r*mu)*a;

	map<double,double> dic;
	double i= mu;
	double value=1;
     int Npts = 200;
     
    double step = upper_bound/Npts;
    for(auto i= step;i<=upper_bound;i=i+step)
        {
        value = Normal(i,mu,sigma);
        dic[i]=value;
        }    


	return dic;
}



void error_model::set_probabilities(double a, double r, double mu, double upper_bound) {
    std::map<double,double> dic1;
    dic1=generate_matrix(a,r,mu,upper_bound);
    std::vector<double> err;
    std::map<double,double> dic= reNormalize(dic1);
    set_deviations(dic);

    for(auto i = dic.begin(); i!=dic.end();++i)
    {   err.push_back(i->second);

    }

      if (_error_dists.empty())
        _error_dists.push_back(err);

    if (_error_dists.size() <= mu)
    {
        _error_dists.resize(200, _error_dists.back());
    }

    _error_dists[mu]=err;
    cout<<endl;
    cout<<mu<<endl;

    cout<<endl;
    cout<<_error_dists.size()<<endl;
    

}

 std::vector<double> error_model::get_probs(double mu,double ux,double scale) const {
  if (mu >= _error_dists.size() && mu <= _max_family_size)
        return _error_dists.back();
    vector<double> retur;
     for(auto i = _error_dists.begin(); i!=_error_dists.end();++i)
    {
        vector<double> x= *i;
        //cout<<mu<<" "<<(x[mu+1]*ux+x[mu]*(1-ux))*scale<<endl;
        retur.push_back((x[mu+1]*ux+x[mu]*(1-ux))*scale);
    }

    return  retur;
}


std::vector<double> error_model::get_epsilons() const {
    set<double> unique_values;
    for (auto& vec : _error_dists)
        unique_values.insert(vec.back());

    vector<double> result(unique_values.size());
    copy(unique_values.begin(), unique_values.end(), result.begin());
    return result;
}

// simple case where we have a single epsilon value in the tree
void error_model::update_single_epsilon(double new_epsilon)
{
    auto epsilons = get_epsilons();
    assert(epsilons.size() == 1);
    map<double, double> replacements;
    replacements[epsilons[0]] = new_epsilon;
    //replace_epsilons(&replacements);
}



