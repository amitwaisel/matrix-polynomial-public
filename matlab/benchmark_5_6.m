function [ all_performance, matlab_performance ] = benchmark_5_6( params )
    params = initialize_params(params);
    all_performance = cell(0);
    matlab_performance = cell(0);

	if ~isfield(params, 'degree_max')
        params.size = 2000;
    end
    if ~isfield(params, 'degree_jump')
        params.size_jump = 250;
    end
    if ~isfield(params, 'degree_from')
        params.size_from = 250;
    end

    params.blocks_algorithm = 'parallel';
    params.schur = 'default';
    params.sylvester = 'recursive';
    params.delta = 0.005;
    params.eigenvalues_sort = 'smart';
    a = -100;
    b = 100;
    
    fprintf('Generating random input matrix with size %d', params.size);
    A = build_unit_eigenvalue_matrix(params.size, params.nworkers);
    fprintf(' ...\n');

    index = 1;
    for co_deg = params.degree_from:params.degree_jump:params.degree_max
    	if isfield(params, 'coefficients')
    		params = rmfield(params, 'coefficients');
    	end
    	params.coefficients = (b-a) .* rand(1,co_deg) + a;
    	coefficients = params.coefficients;

        fprintf('==> Benchmarking poly degree %d for matrix size %d (%s)\n', co_deg, params.size, datestr(datetime('now')));

        all_performance{index, 1} = co_deg;
        all_performance{index, 2} = benchmark_poly_degree(A, params);
        matlab_performance{index, 1} = co_deg;
        matlab_performance{index, 2} = benchmark_matlab_poly_degree(A, params);

        index = index + 1;
	end
    
    fprintf('Done\n');
    detailed_degree_plot(params, all_performance, matlab_performance);
end

function [performance] = benchmark_poly_degree( A, params )
	n = params.size;
	iterations = params.iterations_to_average;
	coefficients = mat2cell(params.coefficients,1,ones(1,numel(params.coefficients)));
	performance = cell(1,5);

    total_time = 0;
    detailed_times = zeros(10,3); % 19
    props = zeros(5,1);
    poly_degree = numel(params.coefficients);
    fprintf('Benchmarking degree %d and size %d (%s): ', poly_degree, n, datestr(datetime('now')));
    fprintf('... ');
    tic;
    for i=1:iterations
        try
            [res, current_time, current_detailed_times, current_props] = MatrixPolynomial(A, coefficients{:}, params);
            total_time = total_time + current_time;
            detailed_times = detailed_times + current_detailed_times;
            props = props + current_props;
        catch err
            fprintf('Error occured for matrix with size %d on iteration %d: %s', n, i, getReport(err));
            rethrow(err)
        end
    end
    performance{1,1} = poly_degree;
    performance{1,2} = total_time/iterations;
    performance{1,3} = detailed_times/iterations;
    performance{1,4} = toc;
    performance{1,5} = props/iterations;
    fprintf('%d (%d)\n', performance{1,2}, performance{1,4});
end

function [performance] = benchmark_matlab_poly_degree( A, params )	
    iterations = params.iterations_to_average;
    performance = cell(1,1);
	n = params.size;

    total_matlab_time = 0;
    fprintf('Benchmarking matlab with polynomial degree %d and size %d (%s): ', numel(params.coefficients), n, datestr(datetime('now')));
    fprintf('... ');
    for i=1:iterations
        tic;
        matlab_polynomial(params.coefficients, A);
        total_matlab_time = total_matlab_time + toc;
    end
    performance{1,1} = total_matlab_time/iterations;
    fprintf('%d\n', performance{1,1});
end

