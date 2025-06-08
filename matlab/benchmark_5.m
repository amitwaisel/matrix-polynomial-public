function [ performance, matlab_performance ] = benchmark_5( params, matlab_params )
    params = initialize_params(params);
    performance = cell(8,1);
        
    if ~isfield(matlab_params, 'matlab_performance')
        [matlab_params.matlab_performance, ~] = benchmark_matlab(params);
    end
    matlab_performance = matlab_params.matlab_performance;
    
    [performance{1,1}, ~] = benchmark_5_0(params, matlab_params);
    [performance{2,1}] = benchmark_5_1(params);
    [performance{3,1}] = benchmark_5_2(params);
%     [performance{4,1}] = benchmark_5_3(params);
    [performance{5,1}] = benchmark_5_4(params);
%     [performance{6,1}, ~] = benchmark_5_5(params);
    [performance{7,1}] = benchmark_5_6(params);
    [performance{8,1}] = benchmark_ps(params);
    
    fprintf('Done\n');
end

