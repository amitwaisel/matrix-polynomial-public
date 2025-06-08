function [fig1,fig2,fig3] = detailed_plot(performance, test_name, matlab_performance, options)
    n = size(performance, 1);
    X = [];
    Y = [];
    
    details_to_show = [];
    if ~exist('options','var')
        options = struct();
    end
        
    for i=1:n
        detailed = cell2mat(performance(i,3));
        if isempty(detailed)
            continue
        end
        ndetails = size(detailed,1);
        details_to_show = [1:(ndetails-4), ndetails-1:ndetails]; % Skip parlett + sylvester details (use total)
        if isfield(options, 'separate_parlett_sylvester') && options.separate_parlett_sylvester
            details_to_show = [1:ndetails-2, ndetails]; % Use parlett_sylvester details, skip total
        end
        if isfield(options, 'optional_details')
            details_to_show = options.optional_details;
        end
        index = performance(i,1);
        d = index{1};
        if isfield(options, 'show_cpu') && options.show_cpu
            detailed = detailed(:,3);
        else
            detailed = detailed(:,1);
        end
        detailed = detailed(details_to_show);
        detailed = (detailed/1e+9);% / (d^2); % Normalize by n
        X = [X, d];
        Y = [Y, detailed];
    end
    
    graph_legend = {'Schur decomposition'       ...
		,'Eigenvalue clustering'       ...
		,'Eigenvalue permutation'      ...
		,'Eigenvalue sorting'          ...
		,'Schur-form reordering'       ...
		,'Polynomial evaluation'       ...
		,'Block Parlett Recurrence'     ...
		,'Sylvester solver'             ...
        ,'Parlett recurrence'           ...
		,'Final multiplication'         ...
        };
    
    graph_legend = graph_legend(details_to_show);
    
    fig1 = figure('Name', test_name, 'Position', [10,10,1280,720]);
    area(X, Y');
    hold on
    
    show_matlab = false;
    if exist('matlab_performance','var') && ~isempty(matlab_performance)
        matlab_performance = cell2mat(matlab_performance(:,1));
        matlab_performance = matlab_performance';% ./ X;
        plot(X, matlab_performance, 'r', 'LineWidth', 3);
        graph_legend{end+1} = 'Naive Horner';
        show_matlab = true;
    end
%     hold on
%     matlab_eval_performance = cell2mat(performance(:,6));
%     plot(X, matlab_eval_performance');
    hold off 
    title_name = get_result_title(test_name);
    xlabel('dimension');
    ylabel('time (s)');
    set_figure_props(title_name, graph_legend);
    
    new_test_name = sprintf('%s', test_name);
    fig2 = figure('Name', new_test_name, 'Position', [10,10,1280,720]);
    plot(X, Y', 'LineWidth', 2);
    hold on;
    new_legend = graph_legend;
    if show_matlab
        plot(X, matlab_performance, 'r', 'LineWidth', 3);
    end
    plot(X, sum(Y,1), 'b', 'LineWidth', 3);
    title_name = get_result_title(new_test_name);
    xlabel('dimension');
    ylabel('time (s)');
    new_legend{end+1} = 'Total';
    if isfield(options, 'ps')
        ps = options.ps;
        ps = cell2mat(ps(X));
        plot(X, ps, 'g', 'LineWidth', 3);
        new_legend{end+1} = 'Paterson-Stockmeyer';
    end
    hold off;
    set_figure_props(title_name, new_legend);
    
    new_test_name = sprintf('%s', test_name);
    fig3 = figure('Name', new_test_name, 'Position', [10,10,1280,720]);
    normY = bsxfun(@times,Y,1./sum(Y,1));
    area(X, normY');
    ylim([0,1]);
    hold on;
    new_legend = graph_legend;
    if show_matlab
        new_legend = new_legend(1:end-1);
    end
    title_name = get_result_title(new_test_name);
    xlabel('dimension');
    ylabel('fraction');
    hold off;
    set_figure_props(title_name, new_legend, true);
end