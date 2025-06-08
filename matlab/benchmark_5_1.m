function [ all_performance ] = benchmark_5_1( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    
    block = 'parallel';
    schur = 'default';
    sort = 'smart';
    sylvester = 'recursive';
    orig_delta = get_best_delta(params.nclusters);
    for delta=[0.1*orig_delta, orig_delta, 5*orig_delta, 100*orig_delta]
        temp_result = benchmark_params(block, sylvester, schur, sort, delta, coefficients, params.nworkers, params);
        nclusters = temp_result(end,5);
        nclusters = nclusters{1}(1);
        result_name = get_result_name(sprintf('x%d eigenvalue clusters', nclusters));
        fprintf('==> Benchmarked %s\n', get_result_title(result_name));
        all_performance.(result_name) = temp_result;
    end
    
    fprintf('Done\n');
    blocks_options = struct(); blocks_options.optional_details = [2, 3, 6, 9, 10];
    plot_performance(all_performance, 'Eigenvalues clustering', [], blocks_options);
end

