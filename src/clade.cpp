#include <iostream>
#include <algorithm>
#include <regex>

#include "doctest.h"
#include "easylogging++.h"

#include "clade.h"
#include "gene_transcript.h"
#include "newick_ape_loader.h"

using namespace std;

clade::clade(const clade& c, clade* parent, std::function<double(const clade& c)> branchlength_setter) {
    _p_parent = parent;
    _taxon_name = c._taxon_name;
    _ape_index = 0;
    if (branchlength_setter)
        _branch_length = branchlength_setter(c);
    else
        _branch_length = c._branch_length;
    _sigma_index = c._sigma_index;
    is_sigma_clade = c.is_sigma_clade;
    _descendants.resize(c._descendants.size());
    transform(c._descendants.begin(), c._descendants.end(), _descendants.begin(), [&](const clade* c) { return new clade(*c, this, branchlength_setter);});

    if (is_root())	// for simplicity, we do not calculate descendantss reverse-level order arrays
        update_reverse_level_order();
}

/* Recursive destructor */
clade::~clade() {

  for (auto i : _descendants) {
    delete i; // recursively delete all descendants
  }
}

const clade *clade::get_parent() const {

  return _p_parent; // memory address
}

double clade::get_branch_length() const
{
    if (is_sigma_clade)
        throw std::runtime_error("Requested branch length from sigma tree");

    return _branch_length;
}

int clade::get_sigma_index() const
{
    if (!is_sigma_clade)
        throw std::runtime_error("Requested sigma index from branch length tree");

    return _sigma_index;
}

/* Adds descendant to vector of descendants */
void clade::add_descendant(clade *p_descendant) {
    
  _descendants.push_back(p_descendant);
  _name_interior_clade();
  if (!is_root()) {
    _p_parent->_name_interior_clade();
  }

  update_reverse_level_order();
}

void clade::update_reverse_level_order()
{
    _reverse_level_order.clear();
    std::stack<const clade*> stack;
    std::queue<const clade*> q;

    q.push(this);
    while (!q.empty())
    {
        /* Dequeue node and make it current */
        auto current = q.front();
        q.pop();
        stack.push(current);

        for (auto i : current->_descendants)
        {
            /* Enqueue child */
            q.push(i);
        }
    }

    while (!stack.empty())
    {
        auto current = stack.top();
        stack.pop();
        _reverse_level_order.push_back(current);
    }

    if (!is_root()) {
        _p_parent->update_reverse_level_order();
    }
}

/* Recursively fills vector of names provided as argument */
void clade::add_leaf_names(vector<string> &vector_names) {

    if (_descendants.empty()) {
        vector_names.push_back(_taxon_name); // base case (leaf), and it starts returning
    }
    else {
        for (auto desc : _descendants) {
            desc->add_leaf_names(vector_names);
        }
    }
}

/* Recursively finds internal nodes, and returns vector of clades */
vector<const clade*> clade::find_internal_nodes() const {

  vector<const clade*> internal_nodes;

  /* Base case: returns empty vector */
  if (is_leaf()) { return internal_nodes; }

  else {
    internal_nodes.push_back(this);

    for (auto i : _descendants) {
      auto descendant = i->find_internal_nodes(); // recursion
      if (!descendant.empty()) { internal_nodes.insert(internal_nodes.end(), descendant.begin(), descendant.end()); }
    }

    return internal_nodes;
  }
}


  /* Recursively find pointer to clade with provided taxon name */
const clade *clade::find_descendant(string some_taxon_name) const {

    const clade *p_descendant = nullptr;
    auto descendant_finder = [some_taxon_name, &p_descendant](const clade *clade) {
        if (clade->get_taxon_name() == some_taxon_name)
            p_descendant = clade;
    };

  apply_prefix_order(descendant_finder);

  return p_descendant;
}
  
/* Finds branch length of clade with provided taxon name. Does so by calling find_descendant(), which recursively searches the tree */
double clade::find_branch_length(string some_taxon_name) {

  auto clade = find_descendant(some_taxon_name);
  if (clade == NULL || clade->is_root()) { return 0; } // guarding against root query

  LOG(INFO) << "Found matching clade";
  return clade->_branch_length;
}

/* Names interior clades, starting from clade of first method call and going up the tree until root */
void clade::_name_interior_clade() {

    vector<string> descendant_names; // vector of names
    add_leaf_names(descendant_names); // fills vector of names
    sort(descendant_names.begin(), descendant_names.end()); // sorts alphabetically (from std)
    _taxon_name.clear(); // resets whatever taxon_name was
    for (auto name : descendant_names) {
        _taxon_name += name;
    }

    if (_p_parent)
        _p_parent->_name_interior_clade();
}

bool clade::is_leaf() const {

  return _descendants.empty();
}

bool clade::is_root() const {

  return get_parent() == NULL;
}

/* Function print_clade_name() is used in conjunction with apply_reverse_level_order and apply_prefix order for debugging purposes.
   e.g.,
   cout << "Tree " << p_tree->get_taxon_name() << " in reverse order: " << endl;
   p_tree->apply_reverse_level_order(print_name)   
*/
void print_clade_name(clade *clade) {
  LOG(INFO) << clade->get_taxon_name() << " (length of subtending branch: " << clade->get_branch_length() << ")" << endl;
}

std::map<std::string, int> clade::get_sigma_index_map()
{
    std::map<std::string, int> node_name_to_sigma_index;
    
    auto fn = [&node_name_to_sigma_index](const clade *c) {
        node_name_to_sigma_index[c->get_taxon_name()] = c->get_sigma_index() - 1; // -1 is to adjust the index offset
    };

    apply_prefix_order(fn);
    return node_name_to_sigma_index;
}

void clade::write_newick(ostream& ost, std::function<void(std::ostream& ost, const clade *c)> write_clade) const
{
    if (is_leaf()) {
        write_clade(ost, this);
    }
    else {
        ost << '(';

        // some nonsense to supress trailing comma
        for (size_t i = 0; i< _descendants.size() - 1; i++) {
            _descendants[i]->write_newick(ost, write_clade);
            ost << ',';
        }

        _descendants[_descendants.size() - 1]->write_newick(ost, write_clade);
        ost << ')';
        write_clade(ost, this);
    }
}

double clade::distance_from_root_to_tip() const
{
    vector<double> candidates;
    apply_prefix_order([&candidates, this](const clade* c) {
        if (c->is_leaf())
        {
            auto p = c;
            double dist = 0;
            while (p)
            {
                dist += p->get_branch_length();
                p = p->get_parent();
            }
            candidates.push_back(dist);
        }
        });

    return *max_element(candidates.begin(), candidates.end());
}

std::ostream& operator<<(std::ostream& ost, const clade& c) {
    if (c.is_leaf())
        ost << c.get_taxon_name();
    ost << "<" << c.get_ape_index() << ">";
    return ost;
}

std::set<double> clade::get_branch_lengths() const
{
    set<double> result;
    auto branch_length_func = [&result](const clade* c) { 
        if (c->get_branch_length() > 0.0)
            result.insert(c->get_branch_length()); 
    };
    apply_prefix_order(branch_length_func);
    return result;
}

void clade::validate_sigma_tree(const clade* p_sigma_tree) const
{
    auto g = [](set<string>& s, const clade* c) {
        s.insert(c->get_taxon_name());
    };
    set<string> my_taxa;
    apply_prefix_order([g, &my_taxa](const clade* c) { g(my_taxa, c);  });

    set<string> sigma_taxa;
    p_sigma_tree->apply_prefix_order([g, &sigma_taxa](const clade* c) { g(sigma_taxa, c);  });

    if (my_taxa != sigma_taxa)
    {
        throw std::runtime_error("The sigma tree structure does not match that of the tree");
    }
}

void clade::apply_to_descendants(const cladefunc& f) const {

    // apply f to direct descendants
    // could replace with apply_prefix_order for functions f that recur through descendants
    //for_each(_descendants.begin(), _descendants.end(), f); // for_each from std
    // for_each apparently passes by value
    for (auto desc : _descendants)
        f(desc);
}

//! apply the functor f to this clade and also to all descendants.
void clade::apply_prefix_order(const cladefunc& f) const {
    std::stack<const clade*> stack;
    stack.push(this);
    while (!stack.empty())
    {
        auto c = stack.top();
        stack.pop();

        // Moving from right to left in the tree because that's what CAFE does
        auto it = c->_descendants.rbegin();
        for (; it != c->_descendants.rend(); ++it)
        {
            stack.push(*it);
        }
        f(c);
    }
}

void clade::normalize()
{
    double sum = 0.0;
    int count = 0;
    for_each(reverse_level_begin(), reverse_level_end(), [&](const clade* node) {
        sum += node->get_branch_length();
        ++count;
    });
    double mean_weight = (count > 0) ? sum / count : 0.0;

    auto recursive_normalize = [&](const clade* c, double mean_weight, auto&& recursive_ref) -> void {
        if (c->is_leaf())
            return;
        for (auto d : c->_descendants) {
            d->_branch_length = d->_branch_length / mean_weight;
            recursive_ref(d, mean_weight, recursive_ref);
        }
    };
    recursive_normalize(this, mean_weight, recursive_normalize);
    LOG(DEBUG) << "Mean branch length in weight tree: " << mean_weight;
}

const clade* common_ancestor(const clade* c1, const clade *c2)
{
    cladevector path;
    while (c1)
    {
        path.push_back(c1);
        c1 = c1->get_parent();
    }
    while (c2)
    {
        if (find(path.begin(), path.end(), c2) != path.end())
            return c2;
        c2 = c2->get_parent();
    }

    return nullptr;
}

double distance(const clade* c1, const clade* c2)
{
    double dist = 0;
    auto p = common_ancestor(c1, c2);
    while (c1 != p)
    {
        dist += c1->get_branch_length();
        c1 = c1->get_parent();
    }
    while (c2 != p)
    {
        dist += c2->get_branch_length();
        c2 = c2->get_parent();
    }
    return dist;
}

clade* parse_newick(std::string newick_string, bool parse_to_sigmas) {

    std::regex tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;");

    auto new_clade = [](clade* p_parent) {
        clade* p_new_clade = new clade();
        if (p_parent != NULL) {
            p_new_clade->_p_parent = p_parent;
        }

        return p_new_clade;
    };

    int lp_count(0), rp_count(0);
    sregex_iterator regex_it(newick_string.begin(), newick_string.end(), tokenizer);
    sregex_iterator regex_it_end;
    clade* p_root_clade = new_clade(NULL);
    p_root_clade->_source_newick = newick_string;

    p_root_clade->is_sigma_clade = parse_to_sigmas; // if user does not provide sigma for root, we need to make the root specifically a sigma clade if we are parsing to sigmas

    clade* p_current_clade = p_root_clade; // current_clade starts as the root

    // The first element below is empty b/c I initialized it in the class body
    for (; regex_it != regex_it_end; regex_it++) {
        /* Checking all regex */
        // cout << regex_it->str() << endl;

        /* Start new clade */
        if (regex_it->str() == "(") {
            /* Checking '(' regex */
            // cout << "Found (: " << regex_it->str() << endl;

            p_current_clade = new_clade(p_current_clade); // move down the tree (towards the present)
            p_current_clade->_p_parent->add_descendant(p_current_clade); // can't forget to add the now current clade to its parent's descendants vector
            lp_count++;
        }

        else if (regex_it->str() == ",") {
            /* Checking ',' regex */
            // cout << "Found ,: " << regex_it->str() << endl;

            /* The if block below is for when the newick notation omits the external parentheses, which is legal */
            if (p_current_clade == p_root_clade) {
                cout << "Found root!" << endl;
                p_root_clade = new_clade(NULL);
                p_current_clade->_p_parent = p_root_clade; // note that get_parent() cannot be used here because get_parent() copies the pointer and it would be the copy that would be assigned p_root_clade... and then the copy would just be thrown away
                p_current_clade->_p_parent->add_descendant(p_current_clade);
            }

            /* Start new clade at same level as the current clade */
            p_current_clade = new_clade(p_current_clade->_p_parent); // move to the side of the tree
            p_current_clade->_p_parent->add_descendant(p_current_clade); // adding current clade as descendant of its parent
        }

        /* Finished current clade */
        else if (regex_it->str() == ")") {
            /* checking ')' regex */
            // cout << "Found ): " << regex_it->str() << endl;

            p_current_clade = p_current_clade->_p_parent; // move up the tree (into the past)
            rp_count++;
        }

        /* Finished newick string */
        else if (regex_it->str() == ";") {
            /* Checking ';' regex */
            // cout << "Found ;: " << regex_it->str() << endl;
            break;
        }

        /* Reading branch length */
        else if (regex_it->str()[0] == ':') {
            /* Checking ':' regex */
            // cout << "Found :: " << regex_it->str() << endl;

            if (parse_to_sigmas)
            {
                int ind = strtol(regex_it->str().substr(1).c_str(), nullptr, 0);
                p_current_clade->_sigma_index = ind;
                p_current_clade->is_sigma_clade = true;
            }
            else
            {
                p_current_clade->_branch_length = atof(regex_it->str().substr(1).c_str()); // atof() converts string into float
                p_current_clade->is_sigma_clade = false;
            }
        }

        /* Reading taxon name */
        else {
            /* Checking species name string regex */
            // cout << "Found species name: " << regex_it->str() << endl;

            p_current_clade->_taxon_name = regex_it->str();
            clade* p_parent = p_current_clade->_p_parent;
            /* If this species has a parent, we need to update the parent's name */
            if (p_parent != NULL) {
                p_parent->_name_interior_clade(); // update parent's name, _name_interior_clade() is a void method
            }
        }
    }

    std::function<void(const clade*)> validator;
    if (p_root_clade->is_sigma_clade)
    {
        // since user is not required to set a sigma index for the root, go ahead and assign it to the first sigma
        // so the rest of the code doesn't get confused
        if (p_root_clade->get_sigma_index() == 0)
            p_root_clade->_sigma_index = 1;

        validator = [](const clade* c) {
            if (c->_sigma_index < 1)
                throw std::runtime_error("Invalid sigma index set for " + c->get_taxon_name());
        };
    }
    else
    {
        validator = [](const clade* c) {
            if (!c->is_root() && c->get_branch_length() <= 0)
                throw std::runtime_error("Invalid branch length set for " + c->get_taxon_name());
        };
    }
    for_each(p_root_clade->reverse_level_begin(), p_root_clade->reverse_level_end(), validator);

    auto order = get_ape_order(p_root_clade);
    std::function<void(clade*)> assign_id;
    
    assign_id = [&order, &assign_id](clade* c) {
        c->_ape_index = distance(order.begin(), find(order.begin(), order.end(), c));
        for_each(c->_descendants.begin(), c->_descendants.end(), assign_id);
    };

    assign_id(p_root_clade);
    return p_root_clade;
}

TEST_CASE("distance_from_root_to_tip")
{
    unique_ptr<clade> p_tree(parse_newick("(A:1,B:3):7"));
    CHECK_EQ(10, p_tree->distance_from_root_to_tip());
}

TEST_CASE("Newick tree is recoverable at the root")
{
    string nwk = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parse_newick(nwk));
    CHECK_EQ(nwk, p_tree->get_source_newick());
}

TEST_CASE("Reverse Level Order touches lowest nodes before parent nodes")
{
    string nwk = "((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;";

    unique_ptr<clade> p_tree(parse_newick(nwk));
    string seq;
    for_each(p_tree->reverse_level_begin(), p_tree->reverse_level_end(), [&seq](const clade* c) {
        seq += c->get_taxon_name();
        });
    CHECK_EQ("BAFCDEGHI", seq);
}

TEST_CASE("common_ancestor")
{
    string nwk = "((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;";
    unique_ptr<clade> p_tree(parse_newick(nwk));
    CHECK_EQ("G", common_ancestor(p_tree->find_descendant("A"), p_tree->find_descendant("C"))->get_taxon_name());
}


TEST_CASE("distance")
{
    string nwk = "((E:0.36,D:0.30)H:1.00,(C:0.85,(A:0.59,B:0.35)F:0.42)G:0.45)I;";
    unique_ptr<clade> p_tree(parse_newick(nwk));
    CHECK_EQ( doctest::Approx(1.86), distance(p_tree->find_descendant("A"), p_tree->find_descendant("C")));
}

TEST_CASE("normalize")
{
    string nwk = "((E:1,D:2)H:3,(C:4,(A:5,B:6)F:7)G:8)I;";
    unique_ptr<clade> p_tree(parse_newick(nwk));
    p_tree->normalize();
    ostringstream ost;
    p_tree->write_newick(ost, [](ostream& ost, const clade* c) {
        ost << c->get_taxon_name() << ":" << c->get_branch_length();
    } );
    CHECK_EQ( "((E:0.25,D:0.5)H:0.75,(C:1,(A:1.25,B:1.5)F:1.75)G:2)I:0", ost.str());
}
