#ifndef io_h
#define io_h

#include <vector>

class clade;
class error_model;
class gene_transcript;
class missing_expression_value;

clade *read_tree(std::string tree_file_path, bool lambda_tree);

void read_gene_transcripts(std::istream& input_file, clade *p_tree, std::vector<gene_transcript>& transcripts);

void read_error_model_file(std::istream& error_model_file, error_model *p_error_model);
void write_error_model_file(std::ostream& ost, error_model& errormodel);

std::vector<std::string> tokenize_str(std::string some_string, char some_delim);

void create_directory(std::string& dir);

std::ostream& operator<<(std::ostream& ost, const gene_transcript& family);

std::ostream& operator<<(std::ostream& ost, const clade& c);

template<typename T>
void write_node_ordered(std::ostream& ost, std::string title, const std::vector<const clade*>& order, std::function<T(const clade* c)> f = nullptr)
{
    ost << title;
    for (auto node : order)
    {
        if (node)
        {
            ost << "\t";
            if (f)
            {
                try
                {
                    ost << f(node);
                }
                catch (missing_expression_value&)
                {
                    ost << "N";
                }
            }
            else
                ost << *node;
        }
    }
    ost << std::endl;
}

#endif
