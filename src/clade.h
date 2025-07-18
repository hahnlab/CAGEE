#ifndef clade_h
#define clade_h

#include <map>
#include <stack>
#include <queue>
#include <string>
#include <functional>
#include <set>

/* Forward declaration of newick_parser class, so class clade can see friend */
class newick_parser; // actual declaration in utils.h

class clade;
using cladefunc = std::function<void(const clade*)>;
class fitch_margoliash;

/*! \brief A Clade represents a node in a tree
*
*  In biology, a clade represents a group of organisms believed to have evolved from a common ancestor.
*  The Clade class has a parent clade, and a list of descendant clades. It can be loaded from a file
*  via the @newick_parser class.
*/
class clade {

    friend clade* parse_newick(std::string newick_string, bool parse_sigmas); // allows newick_parser to set parameter values
    friend class fitch_margoliash; // allows fitch_margoliash to access private members

private:
    clade *_p_parent; // needs to be pointer; instance creates infinite loop
    std::string _taxon_name;
    double _branch_length; // or sigma value
    int _sigma_index;
    bool is_sigma_clade;
    int _ape_index;

    std::vector<clade*> _descendants; // same as above

    /* methods */
    void _name_interior_clade();

    std::vector<const clade*> _reverse_level_order;
    void update_reverse_level_order();

    std::string _source_newick;
public:
    typedef std::vector<const clade*>::const_iterator reverse_level_iterator;
    typedef std::vector<clade*>::const_iterator descendant_iterator;

    /* methods */
    clade() : _p_parent(nullptr), _branch_length(0), _sigma_index(0), is_sigma_clade(false), _ape_index(0) {} // basic constructor

    //! constructor giving taxon name and branch length
    clade(std::string taxon_name, double length) : _p_parent(nullptr), 
        _taxon_name(taxon_name), _branch_length(length), _sigma_index(0), 
        is_sigma_clade(false), _ape_index(0) {}

    clade(const clade& c, clade *parent = nullptr, std::function<double(const clade& c)> branchlength_setter = nullptr);

    ~clade(); // destructor

    //! return the parent clade, NULL if there is none
    const clade *get_parent() const;

    //! Add the descendant clade. Used when constructing a tree
    void add_descendant(clade *p_descendant);

    //!
    void add_leaf_names(std::vector<std::string>& vector_names);

    bool is_leaf() const;

    bool is_root() const;

    double get_branch_length() const;

    double distance_from_root_to_tip() const;

    //! In a multiple sigma situation, returns the index of the sigma associated with this particular clade
    int get_sigma_index() const;

    //! returns descendant nodes of this clade that are not leaves
    std::vector<const clade*> find_internal_nodes() const;

    //! returns a descendant clade by the name
    const clade *find_descendant(std::string some_taxon_name) const;

    double find_branch_length(std::string some_taxon_name);

    std::string get_taxon_name() const { return _taxon_name; }

    void write_newick(std::ostream& ost, std::function<void(std::ostream& ost, const clade *c)> write_clade) const;

    std::map<std::string, int> get_sigma_index_map();

    //! Return a unique list of all brnach lengths for this clade and its descendants
    std::set<double> get_branch_lengths() const;

    /// Checks that the list of node names of the sigma tree matches this one
    /// throw an exception if not
    void validate_sigma_tree(const clade* p_sigma_tree) const;

    //! apply the function f to direct descendants. Does not automatically recurse.
    void apply_to_descendants(const cladefunc& f) const;

    //! apply the function f to this clade and also to all descendants.
    void apply_prefix_order(const cladefunc& f) const;

    reverse_level_iterator reverse_level_begin() const {
        return _reverse_level_order.begin();
    }
    reverse_level_iterator reverse_level_end() const {
        return _reverse_level_order.end();
    }

    descendant_iterator descendant_begin() const {
        return _descendants.begin();
    }
    descendant_iterator descendant_end() const {
        return _descendants.end();
    }

    std::string get_source_newick() const {
        return _source_newick;
    }

    int get_ape_index() const {
        return _ape_index;
    }

    // Divide all nodes by the mean branch length    
    void normalize();
};

template<typename T>
using clademap = std::map<const clade *, T>;

using cladevector = std::vector<const clade *>;

clade* parse_newick(std::string newick_string, bool parse_sigmas);
inline clade* parse_newick(std::string newick_string) { return parse_newick(newick_string, false); }

/// @brief Computes the distance between two clades.
/// The distance is defined as the sum of branch lengths from the first clade to the common 
/// ancestor and from the second clade to the common ancestor.
/// @param c1 
/// @param c2 
/// @return distance between c1 and c2
double distance(const clade* c1, const clade* c2);

/// @brief Finds the common ancestor of two clades.
/// @param c1 
/// @param c2
/// @return pointer to the common ancestor clade, or nullptr if no common ancestor is found
/// Note: This function traverses the parent chain of both clades to find the common ancestor
const clade* common_ancestor(const clade* c1, const clade *c2);

#endif
