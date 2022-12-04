#include <numeric>
#include <set>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "error_model.h"

using namespace std;

error_model::error_model()
{
    _deviations = {};
}

void error_model::set_max_family_size(size_t max_cnt) {
    _max_family_size = max_cnt;
}

void error_model::set_deviations(std::map<int,double> deviations) {
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

double error_model::Normal(int x, int mu, double sigma)
{
double deno,numo,density;
deno = pow((2*M_PI*sigma*sigma),-0.5);
numo = exp(-0.5*pow(((x-mu)/sigma),2));
density = numo*deno;
return density;

}

double error_model::Sum(std::map<int,double> dic)
{
    double sum_=0;
     for(auto i = dic.begin(); i!=dic.end();++i)
{    sum_=sum_+i->second;}

return sum_;
}


std::map<int,double> error_model::reNormalize(std::map<int,double>  dic){
	double sum_values = Sum(dic);
    for(auto i = dic.begin(); i!=dic.end();++i)
{   i->second= i->second/sum_values;
}

	return dic;
}



std::map<int,double> error_model::generate_matrix(double a,double r,double mu,int upper_bound, double ux){
	double sigma= exp(r*mu)*a;
	map<int,double> dic;
	int i=round(mu);
	double value=1;
     int Npts = 200;
	while (round(value*10000.00)/10000 !=0 && i!=upper_bound){
        double nx = (Npts - 1) * (i -0) / double(upper_bound );
        int ix = floor(nx);
		value = Normal(ix,ux*mu,(ux*ux)*sigma);
		dic[i]+=value;
		i=i+1;
    }
	int limit=round(mu)-i;
	i=round(mu)-1;
	value=1;
    while (round(value*100000.00)/100000 !=0 && i!=0){
        double nx = (Npts - 1) * (i -0) / double(upper_bound );
        int ix = floor(nx);
		value = Normal(ix,ux*mu,(ux*ux)*sigma);
		dic[i]+=value;
		i=i-1;
    }

	return dic;
}



void error_model::set_probabilities(double a, double r, size_t mu, double upper_bound,double ux) {
    std::map<int,double> dic1;
    dic1=generate_matrix(a,r,mu,upper_bound,ux);
    std::vector<double> err;
    std::map<int,double> dic= reNormalize(dic1);
    set_deviations(dic);
    for(auto i = dic.begin(); i!=dic.end();++i)
    {   err.push_back(i->second);
    }

    _error_dists[mu]=err;





}

 std::vector<double> error_model::get_probs(size_t mu) const {

    return  _error_dists[mu];
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



