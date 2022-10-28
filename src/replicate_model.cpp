#include "replicate_model.h"
#include "gene_transcript.h"

#include "doctest.h"
#include "easylogging++.h"

void replicate_model::apply() const
{

}

TEST_CASE("replicate_model replaces leaf values by adding together the replicates")
{
	gene_transcript g;
	g.set_expression_value("cow-1", 5.0);
	g.set_expression_value("cow-2", 7.0);
}

