function [ all_performance ] = benchmark_5_4( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    
    block = 'parallel';
    schur = 'default';
    sort = 'smart';
    % delta = 0.005;
    sylvester = 'recursive';
    params.nclusters = 50;
    delta = get_best_delta(params.nclusters);
    for nworkers=[1,2,4,8,14,28,56]
        result_name = get_result_name(sprintf('x%d cores', nworkers));
        fprintf('==> Benchmarking %s\n', result_name);
        params.nworkers = nworkers;
        all_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients, nworkers, params);
    end
   
    fprintf('Done\n');
    plot_performance(all_performance, 'Parallelism over multiple cores');
    plot_cores_performance(all_performance, 'Speedup graph for multiple cores');
    detailed_cpu_plot(all_performance, 'Detailed performance for multiple cores');
end


