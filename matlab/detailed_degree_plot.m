function [] = detailed_degree_plot(params, performance, matlab_performance)
    ntests = size(performance, 1);
    X = [];
    Y = [];
    matlab_times = [];
    d = params.size;
    for i=1:ntests
    	current_degree = performance{i,1};
    	current_performance = performance{i,2};
        detailed = cell2mat(current_performance(1,3));
        ndetails = size(detailed,1);
        details_to_show = [1:(ndetails-4), ndetails-1:ndetails]; % Skip parlett + sylvester details (use total)
        if isempty(detailed)
            continue
        end
        detailed = detailed(:,1);
        detailed = (detailed/1e+9); % Normalize by n
        detailed = detailed(details_to_show);
        X = [X, current_degree];
        Y = [Y, detailed];
        matlab_time = matlab_performance{i, 2};
        matlab_time = matlab_time{1,1};
        matlab_times = [matlab_times, matlab_time];
    end
    test_name = 'Performance by polynomial degree';
    figure('Name', test_name, 'Position', [10,10,1280,720]);
    hold on
    area(X, Y');
    matlab_times = matlab_times';
    plot(X, matlab_times, 'r', 'LineWidth', 3);
    xlim([X(1), X(end)]);
    graph_legend = {'Schur decomposition'       ...
		,'Eigenvalues clustering'       ...
		,'Eigenvalues permutation'      ...
		,'Eigenvalues sorting'          ...
		,'Eigenvlaues reordering'       ...
		,'Polynomial calculation'       ...
...% 		,'Block Parlett Recurrence'     ...
...% 		,'Sylvester solver'             ...
        ,'Parlett Recurrence'           ...
		,'Final multiplication'         ...
        ,'Naive Horner'};
    title_name = get_result_title(test_name);
    xlabel('Polynimial degree');
    ylabel('Calculation time (seconds)');
    hold off 
    set_figure_props(title_name, graph_legend);
    
    new_test_name = sprintf('%s (detailed)', test_name);
    fig2 = figure('Name', new_test_name, 'Position', [10,10,1280,720]);
    hold on;
    plot(X, Y', 'LineWidth', 2);
    new_legend = graph_legend;
    plot(X, matlab_times, 'r', 'LineWidth', 3);
    plot(X, sum(Y,1), 'b', 'LineWidth', 3);
    title_name = get_result_title(new_test_name);
    xlabel('Polynimial degree');
    ylabel('Calculation time (seconds)');
    new_legend{end+1} = 'Total';
    xlim([X(1), X(end)]);
    hold off;
    set_figure_props(title_name, new_legend);
    
    new_test_name = sprintf('%s (normalized)', test_name);
    fig3 = figure('Name', new_test_name, 'Position', [10,10,1280,720]);
    hold on;
    normY = bsxfun(@times,Y,1./sum(Y,1));
    area(X, normY');
    ylim([0,1]);
    xlim([X(1), X(end)]);
    new_legend = graph_legend;
    new_legend = new_legend(1:end-1);
    title_name = get_result_title(new_test_name);
    xlabel('Polynimial degree');
    ylabel('Calculation time (seconds)');
    hold off;
    set_figure_props(title_name, new_legend, true);
end
