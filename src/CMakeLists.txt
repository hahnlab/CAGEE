target_sources(cagee PRIVATE cagee.cpp DiffMat.cpp core.cpp gamma_core.cpp sigma.cpp optimizer_scorer.cpp user_data.cpp)
target_sources(cagee PRIVATE arguments.cpp easylogging++.cpp reconstruction.cpp likelihood_ratio.cpp rootdist_estimator.cpp)
target_sources(cagee PRIVATE base_model.cpp error_model.cpp gene_transcript.cpp matrix_cache.cpp inference_pruner.cpp)
target_sources(cagee PRIVATE execute.cpp gpu_multiplier.cpp blas_multiplier.cpp root_equilibrium_distribution.cpp)
target_sources(cagee PRIVATE clade.cpp gamma.cpp io.cpp optimizer.cpp simulator.cpp newick_ape_loader.cpp)

target_sources(tests PRIVATE cagee.cpp DiffMat.cpp core.cpp gamma_core.cpp sigma.cpp optimizer_scorer.cpp user_data.cpp)
target_sources(tests PRIVATE arguments.cpp easylogging++.cpp reconstruction.cpp likelihood_ratio.cpp rootdist_estimator.cpp)
target_sources(tests PRIVATE base_model.cpp error_model.cpp gene_transcript.cpp matrix_cache.cpp inference_pruner.cpp)
target_sources(tests PRIVATE execute.cpp gpu_multiplier.cpp blas_multiplier.cpp root_equilibrium_distribution.cpp)
target_sources(tests PRIVATE clade.cpp gamma.cpp io.cpp optimizer.cpp simulator.cpp newick_ape_loader.cpp)

target_sources(diffmat_precalc PRIVATE DiffMat.cpp easylogging++.cpp blas_multiplier.cpp)
