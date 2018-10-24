
#include <numeric>
#include <cmath>
#include <getopt.h>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "src/io.h"
#include "src/core.h"
#include "src/gamma_core.h"
#include "src/root_equilibrium_distribution.h"
#include "src/base_model.h"
#include "src/process.h"
#include "src/gene_family_reconstructor.h"
#include "src/matrix_cache.h"
#include "src/gamma_bundle.h"
#include "src/probability.h"
#include "src/family_generator.h"
#include "src/execute.h"
#include "src/user_data.h"
#include "src/optimizer_scorer.h"

TEST_GROUP(GeneFamilies)
{
};

TEST_GROUP(Inference)
{
    void setup()
    {
        srand(10);
    }
};

TEST_GROUP(Simulation)
{
};

TEST_GROUP(Probability)
{
};

TEST_GROUP(Clade)
{
};

TEST_GROUP(Options)
{
    char *values[100];
    int argc;

    void initialize(vector<string> arguments)
    {
        optind = 0;
        argc = arguments.size();
        for (size_t i = 0; i < arguments.size(); ++i)
        {
            values[i] = strdup(arguments[i].c_str());
        }
    }

    void teardown()
    {
        for (size_t i = 0; i < argc; ++i)
        {
            free(values[i]);
        }
    }
};

input_parameters read_arguments(int argc, char *const argv[]);

TEST(Options, input_short)
{
    initialize({ "cafexp", "-ifile" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, input_long)
{
    initialize({ "cafexp", "--infile", "file" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, input_short_space_separated)
{
    initialize({ "cafexp", "-i", "file" });

    auto actual = read_arguments(argc, values);
    STRCMP_EQUAL("file", actual.input_file_path.c_str());
}

TEST(Options, simulate_long)
{
    initialize({ "cafexp", "--simulate=1000" });

    auto actual = read_arguments(argc, values);
    CHECK_EQUAL(1000, actual.nsims);
}

TEST(Options, simulate_short)
{
    initialize({ "cafexp", "-s1000" });

    auto actual = read_arguments(argc, values);
    CHECK_EQUAL(1000, actual.nsims);
}

TEST(Options, cannot_have_space_before_optional_parameter)
{
    try
    {
        initialize({ "cafexp", "-s", "1000" });

        auto actual = read_arguments(argc, values);
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Unrecognized parameter: '1000'", err.what());
    }
}

TEST(Options, must_specify_lambda_and_input_file_for_estimator)
{
    try
    {
        input_parameters params;
        params.fixed_lambda = 0.05;
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Options -l and -i must both be provided an argument.", err.what());
    }
}

TEST(Options, must_not_specify_rootdist_when_simulating)
{
    try
    {
        input_parameters params;
        params.is_simulating = true;
        params.nsims = 5;
        params.rootdist = "foo.txt";
        params.check_input();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Option -s cannot be provided an argument if -f is specified.", err.what());
    }
}

TEST(GeneFamilies, read_gene_families_reads_cafe_files)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);
    LONGS_EQUAL(5, families.at(0).get_species_size("A"));
    LONGS_EQUAL(10, families.at(0).get_species_size("B"));
    LONGS_EQUAL(2, families.at(0).get_species_size("C"));
    LONGS_EQUAL(6, families.at(0).get_species_size("D"));
}

TEST(GeneFamilies, read_gene_families_reads_simulation_files)
{
    std::string str = "#A\n#B\n#AB\n#CD\n#C\n#ABCD\n#D\n35\t36\t35\t35\t36\t34\t34\t1\n98\t96\t97\t98\t98\t98\t98\t1\n";
    std::istringstream ist(str);

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    std::vector<gene_family> families;
    read_gene_families(ist, p_tree, &families);
    LONGS_EQUAL(35, families.at(0).get_species_size("A"));
    LONGS_EQUAL(36, families.at(0).get_species_size("B"));
    LONGS_EQUAL(36, families.at(0).get_species_size("C"));
    LONGS_EQUAL(34, families.at(0).get_species_size("D"));
    delete p_tree;
}

TEST(Inference, infer_processes)
{
    vector<gene_family> families;
    gene_family fam;
    fam.set_species_size("A", 1);
    fam.set_species_size("B", 2);
    families.push_back(fam);
    gene_family fam2;
    fam2.set_species_size("A", 2);
    fam2.set_species_size("B", 1);
    families.push_back(fam2);
    gene_family fam3;
    fam3.set_species_size("A", 3);
    fam3.set_species_size("B", 6);
    families.push_back(fam3);
    gene_family fam4;
    fam4.set_species_size("A", 6);
    fam4.set_species_size("B", 3);
    families.push_back(fam4);

    single_lambda lambda(0.01);
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:1);";
    clade *p_tree = parser.parse_newick();

    base_model core(&lambda, p_tree, &families, 56, 30, NULL, NULL);
    core.start_inference_processes();


    uniform_distribution frq;
    double multi = core.infer_processes(&frq);
    //core.get_likelihoods();
    DOUBLES_EQUAL(41.7504, multi, 0.001);
    delete p_tree;
}

TEST(Inference, uniform_distribution)
{
    std::vector<int> rd(10);
    std::fill(rd.begin(), rd.end(), 1);
    uniform_distribution ef;
    ef.initialize(rd);
    DOUBLES_EQUAL(.1, ef.compute(5), 0.0001);
}

TEST(Inference, gamma_set_alpha)
{
    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
    model.set_alpha(0.5, 3);
}

TEST(Inference, gamma_adjust_family_gamma_membership)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);

    gamma_model model(NULL, NULL, NULL, 0, 5, 0, 0, NULL, NULL);
}

TEST(Inference, gamma)
{
    std::string str = "Desc\tFamily ID\tA\tB\tC\tD\n\t (null)1\t5\t10\t2\t6\n\t (null)2\t5\t10\t2\t6\n\t (null)3\t5\t10\t2\t6\n\t (null)4\t5\t10\t2\t6\n";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);

    single_lambda lambda(0.5);
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();

    std::vector<int> rootdist;
    gamma_model core(&lambda, p_tree, &families, 0, 5, 0, 0, NULL, NULL);
    core.set_max_sizes(148, 122);

    core.start_inference_processes();
    uniform_distribution frq;
    double actual = core.infer_processes(&frq);
    // DOUBLES_EQUAL(56.3469, actual, .0001);

    delete p_tree;
}

TEST(Inference, stash_stream)
{
    family_info_stash stash;
    stash.family_id = "F01";
    stash.lambda_multiplier = 2.5;
    stash.family_likelihood = 3.7;
    stash.posterior_probability = 4.9;

    std::ostringstream ost;
    ost << stash;
    STRCMP_EQUAL("F01\t2.5\t0\t3.7\t4.9\tN/S", ost.str().c_str());

}

TEST(Probability, probability_of_some_values)
{
    matrix_cache calc(0);
    double lambda = 0.05;
    double branch_length = 5;
    DOUBLES_EQUAL(0.0152237, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 5, 9), 0.00001);

    DOUBLES_EQUAL(0.17573, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 9), 0.00001);

    DOUBLES_EQUAL(0.182728, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);

    branch_length = 1;
    DOUBLES_EQUAL(0.465565, calc.get_from_parent_fam_size_to_c(lambda, branch_length, 10, 10), 0.00001);
}

bool operator==(matrix& m1, matrix& m2)
{
    if (m1.size() != m2.size())
        return false;

    for (int i = 0; i < m1.size(); ++i)
    {
        for (int j = 0; j < m1.size(); ++j)
            if (abs(m1.get(i,j) - m2.get(i,j)) > 0.00001)
                return false;
    }

    return true;
}

TEST(Probability, matrices_take_fractional_branch_lengths_into_account)
{
    matrix_cache calc(141);
    single_lambda lambda(0.006335);
    std::set<double> branch_lengths{ 68, 68.7105 };
    calc.precalculate_matrices(get_lambda_values(&lambda), branch_lengths);
    DOUBLES_EQUAL(0.194661, calc.get_matrix(68.7105, 0.006335).get(5,5), 0.00001); // a value 
    DOUBLES_EQUAL(0.195791, calc.get_matrix(68, 0.006335).get(5, 5), 0.00001);
}

TEST(Probability, the_probability_of_going_from_parent_fam_size_to_c)
{
    DOUBLES_EQUAL(0.194661, the_probability_of_going_from_parent_fam_size_to_c(.006335, 68.7105, 5, 5), 0.00001);
}

TEST(Probability, probability_of_matrix)
{
    matrix_cache calc(5);
    single_lambda lambda(0.05);
    std::set<double> branch_lengths{ 5 };
    calc.precalculate_matrices(get_lambda_values(&lambda), branch_lengths);
    matrix actual = calc.get_matrix(5, lambda.get_single_lambda());
    matrix expected(5);
    double values[5][5] = {
    {1,0,0,0,0},
    { 0.2,0.64,0.128,0.0256,0.00512 },
    { 0.04,0.256,0.4608,0.17408,0.0512 },
    { 0.008,0.0768,0.26112,0.36352,0.187392 },
    { 0.0016,0.02048,0.1024,0.249856,0.305562 } };
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            expected.set(i, j, values[i][j]);
    CHECK(actual == expected);

    // a second call should get the same results as the first
    actual = calc.get_matrix(5, lambda.get_single_lambda());
    CHECK(actual == expected);
}

TEST(Probability, get_conditional_distribution_matrix)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lam(0.05);
    matrix_cache cache(15);
    cache.precalculate_matrices(vector<double>{0.05}, set<double>{1});
    auto cd_matrix = get_conditional_distribution_matrix(p_tree.get(), 10, 12, 3, &lam, cache);
    LONGS_EQUAL(10, cd_matrix.size());
}

TEST(Probability, category_selector)
{
    category_selector cs;
    CHECK(cs.empty());
    cs.weighted_cat_draw(5, { 0.2, 0.6 });
    CHECK_FALSE(cs.empty());
    // TODO: Find some way to test the random_device effectively
}

TEST(Inference, base_optimizer_guesses_lambda_only)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> tree(parser.parse_newick());
    single_lambda sl(0.05);

    base_model model(&sl, tree.get(), NULL, 0, 5, NULL, NULL);
    user_data data;
    unique_ptr<optimizer_scorer> opt(model.get_lambda_optimizer(data));
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.565811, guesses[0], 0.0001);
    delete model.get_lambda();
}

TEST(Inference, base_optimizer_guesses_lambda_and_unique_epsilons)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    error_model err;
    err.set_probs(0, { .0, .7, .3 });
    err.set_probs(1, { .4, .2, .4 });
    std::vector<gene_family> fams;
    single_lambda lambda(0.001);
    base_model model(&lambda, p_tree.get(), &fams, 10, 10, NULL, &err);
    user_data data;
    data.p_error_model = &err;
    unique_ptr<optimizer_scorer> opt(model.get_lambda_optimizer(data));
    CHECK(opt);
    CHECK(dynamic_cast<lambda_epsilon_optimizer*>(opt.get()) != NULL);
#ifdef FALSE_MEMORY_LEAK_FIXED
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(3, guesses.size());
    DOUBLES_EQUAL(0.0808301, guesses[0], 0.0001);
    DOUBLES_EQUAL(0.3, guesses[1], 0.0001);
    DOUBLES_EQUAL(0.4, guesses[2], 0.0001);
    vector<double>().swap(guesses);
#endif
    delete model.get_lambda();
}


TEST(Inference, gamma_optimizer_guesses_lambda_and_alpha)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    clade *p_tree = parser.parse_newick();
    single_lambda sl(0.05);

    gamma_model model(&sl, p_tree, NULL, 0, 5, 4, -1, NULL, NULL);
    user_data data;

    unique_ptr<optimizer_scorer> opt(model.get_lambda_optimizer(data));
    CHECK(opt);
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(2, guesses.size());
    DOUBLES_EQUAL(0.218321, guesses[0], 0.0001);
    DOUBLES_EQUAL(0.565811, guesses[1], 0.0001);
    vector<double>().swap(guesses);
   delete p_tree;
   delete model.get_lambda();
}

TEST(Inference, gamma_optimizer_guesses_alpha_if_lambda_provided)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1,(C:1,D:1):1);";
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda sl(0.05);

    gamma_model model(&sl, p_tree.get(), NULL, 0, 5, 4, 0.7, NULL, NULL);
    user_data data;
    data.p_lambda = &sl;
    unique_ptr<optimizer_scorer> opt(model.get_lambda_optimizer(data));
    CHECK(opt);
    auto guesses = opt->initial_guesses();
    LONGS_EQUAL(1, guesses.size());
    DOUBLES_EQUAL(0.565811, guesses[0], 0.0001);
    vector<double>().swap(guesses);
}

TEST(Inference, base_model_reconstruction)
{
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:1):1";
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda sl(0.05);

    std::vector<gene_family> families(1);
    families[0].set_species_size("A", 3);
    families[0].set_species_size("B", 4);

    base_model model(&sl, p_tree.get(), &families, 5, 5, NULL, NULL);

    matrix_cache calc(6);
    calc.precalculate_matrices(get_lambda_values(&sl), set<double>({ 1 }));
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,6 });
    std::unique_ptr<reconstruction> rec(model.reconstruct_ancestral_states(&calc, &dist));

}

TEST(Inference, branch_length_finder)
{
    branch_length_finder finder;
    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    clade *p_tree = parser.parse_newick();
    p_tree->apply_prefix_order(finder);
    LONGS_EQUAL(finder.longest(), 23);
    auto expected = set<double>{ 1, 3, 7, 11, 17, 23 };
    CHECK(finder.result() == expected);
    delete p_tree;

}

TEST(Inference, increase_decrease)
{
    clademap<int> family_size;
    clademap<family_size_change> result;

    newick_parser parser(false);
    parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
    unique_ptr<clade> p_tree(parser.parse_newick());

     auto a = p_tree->find_descendant("A");
     auto b = p_tree->find_descendant("B");
     auto ab = p_tree->find_descendant("AB");
     auto abcd = p_tree->find_descendant("ABCD");

    family_size[a] = 2;
    family_size[b] = 4;
    family_size[ab] = 3;
    family_size[abcd] = 3;
    compute_increase_decrease(family_size, result);

    CHECK(result[a] == Decrease);
    CHECK(result[b] == Increase);
    CHECK(result[ab] == Constant);
}

TEST(Inference, increase_decrease_stream)
{
    ostringstream ost;
    increase_decrease id;
    id.gene_family_id = "1234";
    id.pvalue = 0.02;
    id.change = vector<family_size_change>{ Decrease, Constant, Increase, Increase};
    ost << id;
    STRCMP_EQUAL("1234\t0.02\ty\td\tc\ti\ti\t\n", ost.str().c_str());
}

TEST(Inference, gamma_increase_decrease_stream)
{
    ostringstream ost;
    increase_decrease id;
    id.gene_family_id = "1234";
    id.change = vector<family_size_change>{ Decrease, Constant, Increase, Increase };
    id.category_likelihoods = vector<double>{ .1,.2,.3,.4 };
    id.pvalue = 0.02;
    ost << id;
    STRCMP_EQUAL("1234\t0.02\ty\td\tc\ti\ti\t0.1\t0.2\t0.3\t0.4\t\n", ost.str().c_str());
}

TEST(Inference, precalculate_matrices_calculates_all_lambdas_all_branchlengths)
{
    matrix_cache calc(5);
    std::map<std::string, int> m;
    multiple_lambda lambda(m, vector<double>({ .1, .2, .3, .4 }));
    calc.precalculate_matrices(get_lambda_values(&lambda), set<double>({ 1,2,3 }));
    LONGS_EQUAL(12, calc.get_cache_size());
}

TEST_GROUP(Reconstruction)
{
    gene_family fam;
    unique_ptr<clade> p_tree;
    cladevector order;

    void setup()
    {
        newick_parser parser(false);
        parser.newick_string = "((A:1,B:3):7,(C:11,D:17):23);";
        p_tree.reset(parser.parse_newick());

        fam.set_id("Family5");
        fam.set_species_size("A", 11);
        fam.set_species_size("B", 2);
        fam.set_species_size("C", 5);
        fam.set_species_size("D", 6);

        vector<string> nodes{ "A", "B", "C", "D", "AB", "CD", "ABCD" };
        order.resize(nodes.size());
        const clade *t = p_tree.get();
        transform(nodes.begin(), nodes.end(), order.begin(), [t](string s) { return t->find_descendant(s); });

    }
};

TEST(Reconstruction, gene_family_reconstructor)
{
  single_lambda lambda(0.05);
  fam.set_species_size("Mouse", 3);

  clade leaf("Mouse", 7);

  matrix_cache calc(8);
  calc.precalculate_matrices({ .1 }, set<double>({ 7 }));
  gene_family_reconstructor process(cout, &lambda, 2.0, NULL, 7, 0, &fam, &calc, NULL);
  process(&leaf);
  auto L = process.get_L(&leaf);

  // L holds the probability of the leaf moving from size 3 to size n
  LONGS_EQUAL(8, L.size());
  DOUBLES_EQUAL(0.0, L[0], 0.0001);
  DOUBLES_EQUAL(0.0586679, L[1], 0.0001);
  DOUBLES_EQUAL(0.146916, L[2], 0.0001);
  DOUBLES_EQUAL(0.193072, L[3], 0.0001);
}

TEST(Reconstruction, gene_family_reconstructor_print_reconstruction)
{
    auto a = p_tree->find_descendant("A");
    // cladevector order({ p_tree.get(), a });
    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    gene_family_reconstructor gfc(&fam, p_tree.get(), values);
    ostringstream ost;

    gfc.print_reconstruction(ost, order);
    STRCMP_EQUAL("  TREE Family5 = ((A_11:1,B_2:3)4_8:7,(C_5:11,D_6:17)5_6:23)6_7;\n", ost.str().c_str());
}

TEST(Reconstruction, gamma_bundle_print_reconstruction_prints_value_for_each_category_and_a_summation)
{
    auto a = p_tree->find_descendant("A");
    clademap<int> values;
    values[p_tree.get()] = 7;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[p_tree->find_descendant("AB")] = 8;
    values2[p_tree->find_descendant("CD")] = 6;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };

    ostringstream ost;
    gamma_bundle bundle(v, values2, clademap<family_size_change>(), p_tree.get(), &fam);
    bundle.print_reconstruction(ost, order);
    STRCMP_EQUAL("  TREE Family5 = ((A_11:1,B_2:3)4_0_8:7,(C_5:11,D_6:17)5_0_6:23)6_7_7;\n", ost.str().c_str());

}

TEST(Reconstruction, gamma_bundle_get_increases_decreases)
{
    auto a = p_tree->find_descendant("A");
    auto ab = p_tree->find_descendant("AB");
    order = { p_tree.get(), a, ab };

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[a] = 11;
    values[ab] = 13;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[a] = 11;
    values2[ab] = 13;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };
    clademap<family_size_change> c;
    c[ab] = Decrease;

    ostringstream ost;
    gamma_bundle bundle(v, values2, c, p_tree.get(), &fam);
    auto actual = bundle.get_increases_decreases(order, 0.05);
    DOUBLES_EQUAL(0.05, actual.pvalue, 0.000001);
    LONGS_EQUAL(3, actual.change.size());
    LONGS_EQUAL(Constant, actual.change[0]);
    LONGS_EQUAL(Constant, actual.change[1]);
    LONGS_EQUAL(Decrease, actual.change[2]);
}

TEST(Reconstruction, gamma_model_reconstruction)
{
    auto a = p_tree->find_descendant("A");

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    clademap<double> values2;
    values2[p_tree.get()] = 7;
    values2[p_tree->find_descendant("AB")] = 8;
    values2[p_tree->find_descendant("CD")] = 6;

    auto gfc = new gene_family_reconstructor(&fam, p_tree.get(), values);
    vector<gene_family_reconstructor *> v{ gfc };
    gamma_bundle bundle(v, values2, clademap<family_size_change>(), p_tree.get(), &fam);
    vector<gamma_bundle *> bundles{ &bundle };

    vector<double> multipliers{ 0.13, 1.4 };
    gamma_model_reconstruction gmr(multipliers, bundles);
    std::ostringstream ost;
    gmr.print_reconstructed_states(ost);

    STRCMP_CONTAINS("#NEXUS", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    STRCMP_CONTAINS("  TREE Family5 = ((A_11:1,B_2:3)1_8_8:7,(C_5:11,D_6:17)2_6_6:23)0_7_7;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());

    STRCMP_CONTAINS("BEGIN LAMBDA_MULTIPLIERS;", ost.str().c_str());
    STRCMP_CONTAINS("  0.13;", ost.str().c_str());
    STRCMP_CONTAINS("  1.4;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST(Reconstruction, print_reconstructed_states_empty)
{
    base_model_reconstruction bmr;
    ostringstream ost;
    bmr.print_reconstructed_states(ost);
    STRCMP_EQUAL("", ost.str().c_str());
}

TEST(Reconstruction, print_reconstructed_states_no_print)
{
    auto a = p_tree->find_descendant("A");

    clademap<int> values;
    values[p_tree.get()] = 7;
    values[p_tree->find_descendant("AB")] = 8;
    values[p_tree->find_descendant("CD")] = 6;

    base_model_reconstruction bmr;
    bmr._rec_processes.push_back(new gene_family_reconstructor(&fam, p_tree.get(), values));
    ostringstream ost;
    bmr.print_reconstructed_states(ost);
    STRCMP_CONTAINS("#NEXUS", ost.str().c_str());
    STRCMP_CONTAINS("BEGIN TREES;", ost.str().c_str());
    STRCMP_CONTAINS("  TREE Family5 = ((A_11:1,B_2:3)1_8:7,(C_5:11,D_6:17)2_6:23)0_7;", ost.str().c_str());
    STRCMP_CONTAINS("END;", ost.str().c_str());
}

TEST(Reconstruction, reconstruction_process_internal_node)
{
    single_lambda s_lambda(0.05);
    double multiplier = 2.0;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);

    unique_ptr<lambda> m(s_lambda.multiply(multiplier));
    matrix_cache calc(25);
    calc.precalculate_matrices(get_lambda_values(m.get()), set<double>({ 1, 3, 7, 11, 17, 23 }));
    gene_family_reconstructor process(cout, &s_lambda, multiplier, NULL, 24, 24, &fam, &calc, NULL);

    process(p_tree->find_descendant("A"));
    process(p_tree->find_descendant("B"));

    auto internal_node = p_tree->find_descendant("AB");
    process(internal_node);
    auto L = process.get_L(internal_node);

    // L holds the probability of the leaf moving from size 3 to size n
    LONGS_EQUAL(25, L.size());
    DOUBLES_EQUAL(0.0, L[0], 0.0001);
    DOUBLES_EQUAL(0.00101688, L[1], 0.0001);
    DOUBLES_EQUAL(0.00254648, L[2], 0.0001);
    DOUBLES_EQUAL(0.0033465, L[3], 0.0001);
}

TEST(Inference, gamma_bundle_prune)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda lambda(0.005);
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8, {2,2,2});
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 }, p_tree.get(), &fam);
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,4,3,2,1 });
    matrix_cache cache(11);
    multiple_lambda ml(map<string, int>(), {0.0005, 0.0025});
    cache.precalculate_matrices(get_lambda_values(&ml), set<double>{1, 3, 7});
    CHECK(bundle.prune({ 0.01, 0.05 }, &dist, cache)); 
    auto cat_likelihoods = bundle.get_category_likelihoods();

    LONGS_EQUAL(2, cat_likelihoods.size());
    DOUBLES_EQUAL(-23.3728, log(cat_likelihoods[0]), 0.0001);
    DOUBLES_EQUAL(-17.0086, log(cat_likelihoods[1]), 0.0001);

}

TEST(Inference, gamma_bundle_prune_returns_false_if_saturated)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());
    single_lambda lambda(0.9);
    inference_process_factory factory(cout, &lambda, p_tree.get(), 10, 8, { 2,2,2 });
    factory.set_gene_family(&fam);
    gamma_bundle bundle(factory, { 0.1, 0.5 }, p_tree.get(), &fam);
    uniform_distribution dist;
    dist.initialize({ 1,2,3,4,5,4,3,2,1 });
    matrix_cache cache(11);
    cache.precalculate_matrices({ 0.09, 0.45 }, set<double>{1, 3, 7});

    CHECK(!bundle.prune({ 0.01, 0.05 }, &dist, cache));
}

TEST(Inference, matrix_cache_key_handles_floating_point_imprecision)
{
    set<matrix_cache_key> keys;
    double t = 0.0;
    for (int i = 0; i < 31; i++)
    {
        t += 0.1;
        matrix_cache_key key(1, t, 0.3);
        keys.insert(key);
    }
    LONGS_EQUAL(31, keys.size());

    matrix_cache_key key(1, 3.0, 0.3);
    LONGS_EQUAL(1, keys.count(key));
}

TEST(Inference, birthdeath_rate_with_log_alpha)
{
    // alpha and coeff are derived values from lambda and t
    // (alpha = lambda*t / 1 + lambda*t, coeff = 1 - 2 * alpha);
    DOUBLES_EQUAL(-1.55455, log(birthdeath_rate_with_log_alpha(46, 45, -3.672556, 0.949177)), 0.00001);
    DOUBLES_EQUAL(-2.20436, log(birthdeath_rate_with_log_alpha(44, 46, -2.617970, 0.854098)), 0.00001);
    DOUBLES_EQUAL(-2.39974, log(birthdeath_rate_with_log_alpha(43, 43, -1.686354, 0.629613)), 0.00001);
    DOUBLES_EQUAL(-2.44301, log(birthdeath_rate_with_log_alpha(43, 44, -1.686354, 0.629613)), 0.00001);
    DOUBLES_EQUAL(-1.58253, log(birthdeath_rate_with_log_alpha(13, 14, -2.617970, 0.854098)), 0.00001);

    DOUBLES_EQUAL(0.107, birthdeath_rate_with_log_alpha(40, 42, -1.37, 0.5), .001);
    DOUBLES_EQUAL(0.006, birthdeath_rate_with_log_alpha(41, 34, -1.262, 0.4), .001);

    DOUBLES_EQUAL(0.194661, birthdeath_rate_with_log_alpha(5, 5, -1.1931291703283662, 0.39345841643135504), 0.0001);
}

TEST(Inference, create_one_model_if_lambda_is_null)
{
    input_parameters params;
    params.input_file_path = "foo";
    user_data data;
    auto models = build_models(params, data);
    LONGS_EQUAL(1, models.size());
    CHECK(dynamic_cast<base_model *>(models[0]));
    for (auto m : models)
        delete m;
}

void build_matrix(matrix& m)
{
    m.set(0, 0, 1);
    m.set(0, 1, 2);
    m.set(0, 2, 3);
    m.set(1, 0, 4);
    m.set(1, 1, 5);
    m.set(1, 2, 6);
    m.set(2, 0, 7);
    m.set(2, 1, 8);
    m.set(2, 2, 9);
}

TEST(Probability, matrix_multiply)
{
    matrix m1(3);
    build_matrix(m1);
    vector<double> m2({ 7, 9, 11 });
    auto result = m1.multiply(m2, 0, 2, 0, 2);
    LONGS_EQUAL(3, result.size());

    DOUBLES_EQUAL(58, result[0], .001);
    DOUBLES_EQUAL(139, result[1], .001);
    DOUBLES_EQUAL(220, result[2], .001);

    matrix m3(8);
    m3.set(3, 3, 1);
    m3.set(3, 4, 2);
    m3.set(3, 5, 3);
    m3.set(4, 3, 4);
    m3.set(4, 4, 5);
    m3.set(4, 5, 6);
    m3.set(5, 3, 7);
    m3.set(5, 4, 8);
    m3.set(5, 5, 9);

    result = m3.multiply(m2, 3, 5, 3, 5);
    LONGS_EQUAL(3, result.size());

    DOUBLES_EQUAL(58, result[0], .001);
    DOUBLES_EQUAL(139, result[1], .001);
    DOUBLES_EQUAL(220, result[2], .001);
}

TEST(Probability, error_model_set_probs)
{
    error_model model;
    model.set_probs(0, { .0, .7, .3 });
    model.set_probs(1, { .2, .6, .2 });
    auto vec = model.get_probs(0);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.0, vec[0], 0.00001);
    DOUBLES_EQUAL(0.7, vec[1], 0.00001);
    DOUBLES_EQUAL(0.3, vec[2], 0.00001);

    vec = model.get_probs(1);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);
}

TEST(Probability, error_model_get_epsilon)
{
    error_model model;
    model.set_probs(0, { .0, .7, .3 });
    model.set_probs(1, { .2, .6, .2 });
    model.set_probs(2, { .1, .8, .1 });
    model.set_probs(3, { .2, .6, .2 });

    auto actual = model.get_epsilons();
    LONGS_EQUAL(3, actual.size());
    vector<double> expected{ .1, .2, .3 };
    CHECK(expected == actual);
}

TEST(Probability, error_model_get_epsilon_zero_zero_must_be_zero)
{
    error_model model;
    try
    {
        model.set_probs(0, { 0.4, 0.3, 0.3 });
        CHECK(false);
    }
    catch(runtime_error& err)
    {
        STRCMP_EQUAL("Cannot have a non-zero probability for family size 0 for negative deviation", err.what());
    }
}

TEST(Probability, error_model_rows_must_add_to_one)
{
    error_model model;
    try
    {
        model.set_probs(1, { 0.3, 0.3, 0.3 });
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Sum of probabilities must be equal to one", err.what());
    }
}

TEST(Probability, error_model_replace_epsilons)
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    map<double, double> replacements;
    replacements[.2] = .3;
    model.replace_epsilons(&replacements);

    auto actual = model.get_probs(0);
    DOUBLES_EQUAL(.7, actual[1], 0.0001);
    DOUBLES_EQUAL(.3, actual[2], 0.0001);

    actual = model.get_probs(1);
    DOUBLES_EQUAL(.3, actual[0], 0.0001);
    DOUBLES_EQUAL(.4, actual[1], 0.0001);
    DOUBLES_EQUAL(.3, actual[2], 0.0001);
}

TEST(Probability, read_error_model)
{
    string input = "maxcnt: 10\ncntdiff: -1 0 1\n"
        "0 0.0 0.8 0.2\n"
        "1 0.2 0.6 0.2\n"
        "2 0.2 0.6 0.2\n"
        "3 0.2 0.6 0.2\n"
        "5 0.2 0.6 0.2\n";

    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);
    auto vec = model.get_probs(0);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.0, vec[0], 0.00001);
    DOUBLES_EQUAL(0.8, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);

    vec = model.get_probs(1);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);


    vec = model.get_probs(4);
    LONGS_EQUAL(3, vec.size());
    DOUBLES_EQUAL(0.2, vec[0], 0.00001);
    DOUBLES_EQUAL(0.6, vec[1], 0.00001);
    DOUBLES_EQUAL(0.2, vec[2], 0.00001);
}

TEST(Probability, matrix_cache_warns_on_saturation)
{
    matrix_cache m(10);
    m.precalculate_matrices({ 0.05, 0.01 }, { 25 });
    ostringstream ost;
    m.warn_on_saturation(ost);
    STRCMP_EQUAL("WARNING: Saturated branch using lambda 0.05 on branch length 25\n", ost.str().c_str());
}

TEST(Probability, matrix_is_zero_except_00)
{
    matrix_cache c(10);
    c.precalculate_matrices({ 0.05, 0.01 }, { 25 });
    CHECK(c.get_matrix(25, 0.05).is_zero_except_00());
    CHECK_FALSE(c.get_matrix(25, 0.01).is_zero_except_00());
}

TEST(Inference, build_reference_list)
{
    std::string str = "Desc\tFamily ID\tA\tB\n"
        "\t (null)1\t5\t10\n"
        "\t (null)2\t5\t7\n"
        "\t (null)3\t5\t10\n"
        "\t (null)4\t5\t7\n";
    std::istringstream ist(str);
    std::vector<gene_family> families;
    read_gene_families(ist, NULL, &families);
    auto actual = build_reference_list(families);
    vector<int> expected({ 0, 1, 0, 1 });
    LONGS_EQUAL(expected.size(), actual.size());

}

TEST(Inference, prune)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family fam;
    fam.set_species_size("A", 3);
    fam.set_species_size("B", 6);
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);
    inference_process process(ost, &lambda, 1.5, p_tree.get(), 20, 20, &fam, { 1,2,3 }, NULL);
    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });
    auto actual = process.prune(cache);

    vector<double> log_expected{ -17.2771, -10.0323 , -5.0695 , -4.91426 , -5.86062 , -7.75163 , -10.7347 , -14.2334 , -18.0458 , 
        -22.073 , -26.2579 , -30.5639 , -34.9663 , -39.4472 , -43.9935 , -48.595 , -53.2439 , -57.9338 , -62.6597 , -67.4173 };

    LONGS_EQUAL(log_expected.size(), actual.size());
    for (size_t i = 0; i<log_expected.size(); ++i)
    {
        DOUBLES_EQUAL(log_expected[i], log(actual[i]), 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_leaf_nodes_correctly)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    likelihood_computer pruner(20, 20, &lambda, family, cache, NULL);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    auto A = p_tree->find_descendant("A");
    pruner(A);
    auto actual = pruner.get_likelihoods(A);

    vector<double> expected{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i<expected.size(); ++i)
    {
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }

    auto B = p_tree->find_descendant("B");
    pruner(B);
    actual = pruner.get_likelihoods(B);

    expected = { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i<expected.size(); ++i)
    {
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_root_nodes_correctly)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    likelihood_computer pruner(20, 20, &lambda, family, cache, NULL);
    cache.precalculate_matrices({ 0.03 }, { 1.0,3.0,7.0 });

    auto AB = p_tree->find_descendant("AB");
    try
    {
        pruner(AB);
        CHECK(false);   
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Child node probabilities not calculated", err.what());
    }

    pruner(p_tree->find_descendant("A"));
    pruner(p_tree->find_descendant("B"));
    pruner(AB);

    auto actual = pruner.get_likelihoods(AB);

    vector<double> log_expected{ -19.7743, -11.6688, -5.85672, -5.66748, -6.61256, -8.59725, -12.2301, -16.4424, -20.9882, -25.7574, 
        -30.6888, -35.7439, -40.8971, -46.1299, -51.4289, -56.7837, -62.1863, -67.6304, -73.1106, -78.6228
    };

    LONGS_EQUAL(log_expected.size(), actual.size());
    for (size_t i = 0; i<log_expected.size(); ++i)
    {
        DOUBLES_EQUAL(log_expected[i], log(actual[i]), 0.0001);
    }
}

TEST(Inference, likelihood_computer_sets_leaf_nodes_from_error_model_if_provided)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";

    gene_family family;
    family.set_species_size("A", 3);
    family.set_species_size("B", 6);

    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.03);

    matrix_cache cache(21);
    cache.precalculate_matrices({ 0.045 }, { 1.0,3.0,7.0 });

    string input = "maxcnt: 20\ncntdiff: -1 0 1\n"
        "1 0.2 0.6 0.2\n"
        "20 0.2 0.6 0.2\n";
    istringstream ist(input);
    error_model model;
    read_error_model_file(ist, &model);

    likelihood_computer pruner(20, 20, &lambda, family, cache, &model);

    auto A = p_tree->find_descendant("A");
    pruner(A);
    auto actual = pruner.get_likelihoods(A);

    vector<double> expected{ 0, 0, 0.2, 0.6, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    LONGS_EQUAL(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i)
    {
        //cout << actual[i] << endl;
        DOUBLES_EQUAL(expected[i], actual[i], 0.0001);
    }
}

TEST(Clade, get_lambda_index_throws_from_branch_length_tree)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    LONGS_EQUAL(7, p_tree->get_branch_length());
    try
    {
        p_tree->get_lambda_index();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Requested lambda index from branch length tree", err.what());
    }

}

TEST(Clade, get_branch_length_throws_from_lambda_tree)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    parser.parse_to_lambdas = true;
    unique_ptr<clade> p_tree(parser.parse_newick());
    LONGS_EQUAL(7, p_tree->get_lambda_index());

    try
    {
        p_tree->get_branch_length();
        CHECK(false);
    }
    catch (runtime_error& err)
    {
        STRCMP_EQUAL("Requested branch length from lambda tree", err.what());
    }

}

TEST(Inference, multiple_lambda_returns_correct_values)
{
    ostringstream ost;
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    map<string, int> key;
    key["A"] = 5;
    key["B"] = 3;
    multiple_lambda ml(key, { .03, .05, .07, .011, .013, .017 });
    DOUBLES_EQUAL(.017, ml.get_value_for_clade(p_tree->find_descendant("A")), 0.0001);
    DOUBLES_EQUAL(.011, ml.get_value_for_clade(p_tree->find_descendant("B")), 0.0001);
}

class mock_optimizer : public optimizer_scorer
{
    double _initial;
    // Inherited via optimizer
    virtual std::vector<double> initial_guesses() override
    {
        return std::vector<double>{_initial};
    }
    virtual double calculate_score(double * values) override
    {
        return 0.0;
    }
public:
    mock_optimizer(double initial) : _initial(initial)
    {

    }

};

class mock_model : public model {
    // Inherited via model
    virtual simulation_process * create_simulation_process(int family_number) override
    {
        return new simulation_process(cout, _p_lambda, 1.0, _p_tree, 0, 0, vector<int>(), family_number, NULL);
    }
    virtual void start_inference_processes() override
    {
    }
    virtual double infer_processes(root_equilibrium_distribution * prior) override
    {
        return 0.0;
    }
    virtual std::string name() override
    {
        return "mockmodel";
    }
    virtual void write_family_likelihoods(std::ostream & ost) override
    {
    }
    virtual reconstruction* reconstruct_ancestral_states(matrix_cache * p_calc, root_equilibrium_distribution * p_prior) override
    {
    }
    virtual optimizer_scorer * get_lambda_optimizer(user_data& data) override
    {
        return nullptr;
    }
public:
    mock_model() : model(NULL, NULL, NULL, 0, 0, NULL)
    {

    }
    void set_lambda(lambda * lambda)
    {
        _p_lambda = lambda;
    }
    void set_tree(clade * tree)
    {
        _p_tree = tree;
    }

    // Inherited via model
    virtual void prepare_matrices_for_simulation(matrix_cache & cache) override
    {
        branch_length_finder lengths;
        _p_tree->apply_prefix_order(lengths);
        cache.precalculate_matrices(get_lambda_values(_p_lambda), lengths.result());
    }
};

TEST(Simulation, simulate_processes_sets_roots_less_than_100_without_rootdist)
{
    srand(10);
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lam(0.05);
    vector<trial *> results(3);
    mock_model model;
    model.set_lambda(&lam);
    model.set_tree(p_tree.get());
    model.simulate_processes(results);
    LONGS_EQUAL(3, results.size());
    CHECK(results[0]->at(p_tree.get()) < 100);
    CHECK(results[1]->at(p_tree.get()) < 100);
    CHECK(results[2]->at(p_tree.get()) < 100);
    for (auto r : results)
        delete r;
}

TEST(Inference, model_vitals)
{
    mock_model model;
    single_lambda lambda(0.05);
    model.set_lambda(&lambda);
    std::ostringstream ost;
    model.write_vital_statistics(ost, 0.01);
    STRCMP_EQUAL("Model mockmodel Result: 0.01\nLambda:            0.05\n", ost.str().c_str());
}

class increasable {
    clade *_p_tree;
public:
    increasable(clade *p_tree) : _p_tree(p_tree)
    {

    }
    std::vector<const clade *> get_taxa()
    {
        vector<string> taxa{ "A", "B", "AB" };
        vector<const clade*> result(3);
        transform(taxa.begin(), taxa.end(), result.begin(), [this](string taxon)->const clade * {
            return _p_tree->find_descendant(taxon);});
        return result;
    }
    increase_decrease get_increases_decreases(std::vector<const clade *>& order, double pvalue)
    {
        increase_decrease result;
        result.gene_family_id = "myid";
        result.pvalue = pvalue;
        result.change = vector<family_size_change>{ Increase, Decrease, Constant };
        return result;
    }
};


TEST(Inference, print_increases_decreases_by_family)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    vector<increasable*> vec;
    vector<double> pvalues;
    ostringstream empty;
    ::print_increases_decreases_by_family(empty, vec, pvalues);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    unique_ptr<increasable> ic(new increasable(p_tree.get()));
    vec.push_back(ic.get());
    pvalues.push_back(0.07);
    ostringstream ost;
    ::print_increases_decreases_by_family(ost, vec, pvalues);
    STRCMP_CONTAINS("#FamilyID\tpvalue\t*\tA\tB\tAB", ost.str().c_str());
    STRCMP_CONTAINS("myid\t0.07\tn\ti\td\tc", ost.str().c_str());
}

TEST(Inference, print_increases_decreases_by_clade)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    vector<increasable*> vec;
    ostringstream empty;
    ::print_increases_decreases_by_clade(empty, vec);
    STRCMP_CONTAINS("No increases or decreases recorded", empty.str().c_str());

    unique_ptr<increasable> ic(new increasable(p_tree.get()));
    vec.push_back(ic.get());
    ostringstream ost;
    ::print_increases_decreases_by_clade(ost, vec);
    STRCMP_CONTAINS("#Taxon_ID\tIncrease/Decrease", ost.str().c_str());
    STRCMP_CONTAINS("A\t1/0", ost.str().c_str());
    STRCMP_CONTAINS("B\t0/1", ost.str().c_str());
}

TEST(Inference, lambda_epsilon_optimizer)
{
    const double initial_epsilon = 0.01;
    error_model err;
    err.set_probs(0, { .0, .99, initial_epsilon });
    err.set_probs(1, { initial_epsilon, .98, initial_epsilon });

    mock_model model;

    single_lambda lambda(0.05);
    lambda_epsilon_optimizer optimizer(&model, &err, NULL, &lambda, 10);
    optimizer.initial_guesses();
    vector<double> values = { 0.05, 0.06 };
    optimizer.calculate_score(&values[0]);
    auto actual = err.get_probs(0);
    vector<double> expected{ 0, .94, .06 };
    CHECK(expected == actual);

    values[1] = 0.04;
    optimizer.calculate_score(&values[0]);
    actual = err.get_probs(0);
    expected = { 0, .96, .04 };
    CHECK(expected == actual);
}

TEST(Inference, gamma_optimizer)
{
    gamma_optimizer optimizer(NULL, NULL);
    auto initial = optimizer.initial_guesses();
    DOUBLES_EQUAL(0.565811, initial[0], 0.00001);
}

TEST(Inference, gamma_lambd_optimizer)
{
    // TODO: Add families and a tree to the Inference test group, since they are required for pretty much everything
    // Remove p_tree pointer from gamma_lambda_optimizer as it is only used to find the longest branch length
    single_lambda lambda(0.05);

    vector<gene_family> families;
    gene_family fam;
    fam.set_species_size("A", 1);
    fam.set_species_size("B", 2);
    families.push_back(fam);

    uniform_distribution frq;

    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    gamma_model m(&lambda,p_tree.get(), &families, 10, 10, 4, 0.25, NULL, NULL);
    gamma_lambda_optimizer optimizer(p_tree.get(), &lambda, &m, &frq);
    vector<double> values{ 0.01, 0.25 };
    optimizer.calculate_score(&values[0]);
}

TEST(Simulation, simulation_process_max_family_size_is_twice_max_rootdist)
{
    vector<int> rootdist{ 3, 5, 9, 11, 4 };
    simulation_process p(cout, NULL, 1, NULL, 0, 0, rootdist, 0, NULL);
    LONGS_EQUAL(22,p.get_max_family_size_to_simulate());
}

TEST(Simulation, base_prepare_matrices_for_simulation_creates_matrix_for_each_branch)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    base_model b(&lam, p_tree.get(), NULL, 0, 0, NULL, NULL);
    matrix_cache m(25);
    b.prepare_matrices_for_simulation(m);
    LONGS_EQUAL(3, m.get_cache_size());
}

TEST(Simulation, gamma_prepare_matrices_for_simulation_creates_matrix_for_each_branch_and_category)
{
    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    single_lambda lam(0.05);
    unique_ptr<clade> p_tree(parser.parse_newick());
    gamma_model g(&lam, p_tree.get(), NULL, 0, 0, 2, 0.5, NULL, NULL);
    matrix_cache m(25);
    g.prepare_matrices_for_simulation(m);
    LONGS_EQUAL(6, m.get_cache_size());
}

TEST(Simulation, random_familysize_setter_without_error_model)
{
    srand(10);

    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.05);
    trial t;
    matrix_cache cache(12);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    random_familysize_setter setter(&t, 10, &lambda, NULL, cache);
    auto b = p_tree->find_descendant("B");
    setter(b);

    LONGS_EQUAL(0, t[b]);

    t[p_tree.get()] = 5;
    setter(b);
    LONGS_EQUAL(5, t[b]);
}

TEST(Simulation, random_familysize_setter_select_size)
{
    matrix_cache cache(50);
    cache.precalculate_matrices({ 0.05, 0.02 }, { 10, 20 });
    random_familysize_setter setter(NULL, 50, NULL, NULL, cache);
    LONGS_EQUAL(2, setter.select_size(5, 0.05, 10, 0.1));
    LONGS_EQUAL(8, setter.select_size(8, 0.05, 10, 0.5));
    LONGS_EQUAL(15, setter.select_size(8, 0.02, 20, 0.99));
//    CHECK(FALSE);
}

TEST(Simulation, random_familysize_setter_with_error_model)
{
    srand(10);

    const int family_size = 10;
    const double initial_epsilon = 0.2;
    error_model err;
    err.set_probs(0, { .0, .8, initial_epsilon });
    err.set_probs(1, { initial_epsilon, .6, initial_epsilon });
    err.set_probs(family_size, { initial_epsilon, .6, initial_epsilon });

    newick_parser parser(false);
    parser.newick_string = "(A:1,B:3):7";
    unique_ptr<clade> p_tree(parser.parse_newick());

    single_lambda lambda(0.05);
    trial t;
    matrix_cache cache(20);
    cache.precalculate_matrices({ 0.05 }, { 3 });
    random_familysize_setter setter(&t, family_size, &lambda, &err, cache);
    auto b = p_tree->find_descendant("B");
    setter(b);

    LONGS_EQUAL(0, t[b]);

    t[p_tree.get()] = 5;
    setter(b);
    LONGS_EQUAL(4, t[b]);
}

TEST(Simulation, executor)
{
    input_parameters params;
    user_data ud;
    unique_ptr<action> act(get_executor(params, ud));
    CHECK(dynamic_cast<const estimator *>(act.get()))

    params.chisquare_compare = true;
    unique_ptr<action> act2(get_executor(params, ud));
    CHECK(dynamic_cast<const chisquare_compare *>(act2.get()))
}

void init_lgamma_cache();

int main(int ac, char** av)
{
    init_lgamma_cache();
    return CommandLineTestRunner::RunAllTests(ac, av);
}