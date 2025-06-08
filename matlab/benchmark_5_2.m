function [ all_performance ] = benchmark_5_2( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    
    schur = 'default';
    sort = 'smart';
    sylvester = 'recursive';
    delta = 0.005;
    for block_cell = {'parallel', 'parlett', 'sequential'}
        block = block_cell{1};
        result_name = get_result_name(sprintf('%s block algorithm', block));
        fprintf('==> Benchmarking %s\n', result_name);
        all_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients, params.nworkers, params);
    end
    
    fprintf('Done\n');
    parlett_options = struct(); parlett_options.optional_details = [6 9 10];
    plot_performance(all_performance, 'Parlett Recurrence technique', [], parlett_options);
    plot_specific_performance(all_performance, 'Parlett recurrence technique', 9)
end

