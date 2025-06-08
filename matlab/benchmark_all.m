function [ all_performance, matlab_performance ] = benchmark_all( params, matlab_params )
    if ~isfield(params, 'size')
        params.size = 100;
    end
    if ~isfield(params, 'iterations_to_average')
        params.iterations_to_average = 10;
    end
    if ~isfield(params, 'delta')
        params.delta = 0.1;
    end
    if ~isfield(params, 'nworkers')
        params.nworkers = 8;
    end
%     if ~isfield(params, 'coefficients')
%         params.coefficients = [1, -3, 19, -4, 0, 0, 9];
%     end
    short_coefficients = [1, -3, 19, -4, 0, 0, 9];
    a = -100;
    b = 100;
    long_coefficients = (b-a) .* rand(1,800) + a;
    coefficients = struct();
    coefficients.long = long_coefficients;
    coefficients.short = short_coefficients;
    all_performance = struct();
    matlab_performance = struct();
    coefficient_types = fieldnames(coefficients);
    for co_type = 1:length(coefficient_types)
        co = coefficient_types{co_type};
        if ~isfield(matlab_params, 'matlab_performance')
            matlab_params.matlab_performance = struct();
        end
        for block_cell = {'parallel'} %'parlett', 'sequential'}
            block = block_cell{1};
            for schur_cell = {'default'} %{'default', 'hassenberg'}
                schur = schur_cell{1};
                for sort_cell = {'smart'} %{'bubble', 'smart'}
                    sort = sort_cell{1};
                    for sylvester_cell = {'recursive'} %{'intel', 'lapack', 'recursive'}
                        sylvester = sylvester_cell{1};
                        for delta=[0.05] %,0.1]
                            result_name = sprintf('%s_%s_%s_%s_%s__%f', co, block, schur, sort, sylvester, delta);
                            result_name = strrep(result_name, '.', '_');
                            fprintf('==> Benchmarking %s\n', result_name);
                            all_performance.(result_name) = benchmark_params(block, sylvester, schur, sort, delta, coefficients.(co), params);
                            
                            
                            if ~isfield(matlab_params.matlab_performance, result_name)
                                matlab_params.matlab_performance.(result_name) = benchmark_matlab(coefficients.(co), params);
                            end
                            
                            matlab_performance.(result_name) = matlab_params.matlab_performance.(result_name);
                        end
                    end
                end
            end
        end
    end
    fprintf('Done\n');
    plot_performance(all_performance, matlab_performance);
end


