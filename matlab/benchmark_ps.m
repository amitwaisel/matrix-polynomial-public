function [ all_performance ] = benchmark_ps( params )
    params = initialize_params(params);
    coefficients = params.coefficients;
    all_performance = struct();
    for ps_cell = {'vanloan', 'parallel', 'sequential'}
        ps_type = ps_cell{1};
        result_name = get_result_name(sprintf('%s polynomial calculation', ps_type));
        fprintf('==> Benchmarking %s\n', result_name);
        all_performance.(result_name) = benchmark_ps_impl(coefficients, ps_type, params);
    end
    fprintf('Done\n');
    plot_ps_performance(all_performance);
end

function [performance] = benchmark_ps_impl(co, ps_type, params)
    if ~isfield(params, 'size')
        params.size = 100;
    end
    if ~isfield(params, 'size_jump')
        params.size_jump = 1;
    end
    if ~isfield(params, 'size_from')
        params.size_from = 2;
    end
    
    n = params.size;
    iterations = params.iterations_to_average;
    performance = cell(n,3);
    coefficients = mat2cell(co,1,ones(1,numel(co)));
    params.delta = 10;
    for s = params.size_from:params.size_jump:n
        total_time = 0;
        fprintf('Benchmarking size %d: ', s);
        A = build_unit_eigenvalue_matrix(s, params.nclusters);
        fprintf('... ');
        tic;
        for i=1:iterations
            try
                [res, current_time] = PolynomialTest(A, coefficients{:}, params);
                total_time = total_time + current_time;
            catch err
                fprintf('Error occured for matrix with size %d on iteration %d', s, i);
                name = sprintf('C:\\temp\\A_%d_%d.mat', s, s);
                save(name, 'A');
                rethrow(err)
            end
        end
        performance{s,1} = s;
        performance{s,2} = total_time/iterations;
        performance{s,3} = toc;
        fprintf('%d (%d)\n', performance{s,2}, performance{s,3});
    end
end


