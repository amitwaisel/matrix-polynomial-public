%DEPRECATED
function [ all_performance ] = benchmark_5_3( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    
    block = 'parallel';
    schur = 'default';
    sort = 'smart';
    delta = 0.005;
    for sylvester_cell = {'recursive', 'intel', 'lapack'}
        sylvester = sylvester_cell{1};
        result_name = get_result_name(sprintf('%s Sylvester solver', sylvester));
        fprintf('==> Benchmarking %s\n', result_name);
        all_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients, params.nworkers, params);
    end
    
    fprintf('Done\n');
    plot_performance(all_performance, 'Sylvester solver types');
end
