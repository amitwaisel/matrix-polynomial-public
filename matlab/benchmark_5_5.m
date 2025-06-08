%DEPRECATED
function [ all_performance ] = benchmark_5_5( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    
    block = 'parallel';
    schur = 'default';
    sylvester = 'recursive';
    delta = 0.005;
    for sort_cell={'smart','bubble'}
        sort=sort_cell{1};
        result_name = get_result_name(sprintf('%s sort for eigenvalues', sort));
        fprintf('==> Benchmarking %s\n', result_name);
        all_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients, params.nworkers, params);
    end
   
    fprintf('Done\n');
    sort_options = struct(); sort_options.optional_details = [4,5];
    plot_performance(all_performance, 'Eigenvalues sorting', [], sort_options);
end
