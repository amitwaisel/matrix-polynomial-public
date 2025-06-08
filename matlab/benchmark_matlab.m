function [ performance, params ] = benchmark_matlab( params )
    params = initialize_params(params);
    co = params.coefficients;
    n = params.size;
    nclusters = params.nclusters;
    iterations = params.iterations_to_average;
    performance = cell(n,1);
    for s = params.size_from:params.size_jump:n
        total_matlab_time = 0;
        fprintf('Benchmarking matlab with polynomial degree %d and size %d: ', size(co,2), s);
        A = build_unit_eigenvalue_matrix(s, nclusters);
        fprintf('... ');
        for i=1:iterations
            tic;
            matlab_polynomial(co, A);
            total_matlab_time = total_matlab_time + toc;
        end
        performance{s,1} = total_matlab_time/iterations;
        fprintf('%d\n', performance{s,1});
    end
end

