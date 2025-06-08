function [ performance ] = benchmark( params )
    if ~isfield(params, 'size')
        params.size = 100;
    end
    if ~isfield(params, 'size_jump')
        params.size_jump = 1;
    end
    if ~isfield(params, 'size_from')
        params.size_from = 2;
    end
    if ~isfield(params, 'iterations_to_average')
        params.iterations_to_average = 10;
    end
    if ~isfield(params, 'delta')
        params.delta = 0.1;
    end
    if ~isfield(params, 'blocks_algorithm')
        params.blocks_algorithm = 'parlett';
    end
    if ~isfield(params, 'sylvester')
        params.sylvester = 'recursive';
    end
    if ~isfield(params, 'nworkers')
        params.nworkers = 8;
    end
    if ~isfield(params, 'schur')
        params.schur = 'default';
    end
    if ~isfield(params, 'eigenvalues_sort')
        params.eigenvalues_sort = 'eigenvalues_sort';
    end
    if ~isfield(params, 'coefficients')
        params.coefficients = [1, -3, 19, -4, 0, 0, 9];
    end    
    if ~isfield(params, 'nclusters')
        params.nclusters = params.nworkers;
    end
    
    n = params.size;
    nclusters = params.nclusters;
    iterations = params.iterations_to_average;
    coefficients = mat2cell(params.coefficients,1,ones(1,numel(params.coefficients)));
    performance = cell(n,5);
    for s = params.size_from:params.size_jump:n
        total_time = 0;
        detailed_times = zeros(10,3); % 19
        props = zeros(5,1);
        fprintf('Benchmarking size %d with %d clusters and delta %f at %s: ', s, params.nclusters, params.delta, datestr(datetime('now')));
        A = build_unit_eigenvalue_matrix(s, nclusters);
        fprintf('... %s ... ', datestr(datetime('now')));
        tic;
        for i=1:iterations
            try
                [res, current_time, current_detailed_times, current_props] = MatrixPolynomial(A, coefficients{:}, params);
                total_time = total_time + current_time;
                detailed_times = detailed_times + current_detailed_times;
                props = props + current_props;
            catch err
                fprintf('Error occured for matrix with size %d on iteration %d: %s', s, i, getReport(err));
                name = sprintf('C:\\temp\\A_%d_%d.mat', s, s);
                save(name, 'A');
                rethrow(err)
            end
        end
        performance{s,1} = s;
        performance{s,2} = total_time/iterations;
        performance{s,3} = detailed_times/iterations;
        performance{s,4} = toc;
        performance{s,5} = props/iterations;
        num_blocks = performance{s,5};
        num_blocks = num_blocks(1,1);
        fprintf('%d (%d) [%d clusters]\n', performance{s,2}, performance{s,4}, num_blocks);
    end
end


