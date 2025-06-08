function [ all_performance, matlab_performance ] = benchmark_5_0( params, matlab_params )
    all_performance = struct();
    pause on
    block = 'parallel';
    schur = 'default';
    sort = 'smart';
    sylvester = 'recursive';
    delta=0.005;
    for degree=[20, 2000]
        params.degree = degree;
        if isfield(params, 'coefficients')
            params = rmfield(params, 'coefficients');
        end
        params = initialize_params(params);
        coefficients = params.coefficients;
        current_performance = struct();
        
        delta = get_best_delta(params.nclusters);

        result_name = get_result_name(sprintf('degree %d', params.degree));
        fprintf('==> Benchmarking %s\n', result_name);
        current_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients, params.nworkers, params);
        all_performance.(result_name) = current_performance.(result_name);
        fprintf('Done %s at %s, waiting to start Matlab benchmark...', result_name, datestr(datetime('now')));
        pause(5);
        matlab_performance_name = sprintf('matlab_performance_%d', degree);
        if ~isfield(matlab_params, matlab_performance_name)
            matlab_params.(matlab_performance_name) = benchmark_matlab(params);
        end
        matlab_performance.(matlab_performance_name) = matlab_params.(matlab_performance_name);

        fprintf('Done\n');
        plot_performance(current_performance, 'Baseline variant', matlab_performance.(matlab_performance_name));
    end
end
