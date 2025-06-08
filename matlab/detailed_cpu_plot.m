function [] = detailed_cpu_plot(all_performance, test_name)
    tests = fieldnames(all_performance);
    tests_names = cellfun(@get_result_title, tests, 'un', 0);
    ntests = length(tests);
    myColors = jet(ntests);
    
    graph_legend = {'Schur decomposition'       ...
		,'Eigenvalues clustering'       ...
		,'Eigenvalues permutation'      ...
		,'Eigenvalues sorting'          ...
		,'Eigenvalues reordering'       ...
		,'Polynomial calculation'       ...
		,'Block Parlett Recurrence'     ...
		,'Sylvester solver'             ...
        ,'Parlett Recurrence'           ...
		,'Final multiplication'         ...
        };
    
    for detail_index=1:size(graph_legend,2)
        new_test_name = sprintf('%s (%s)', test_name, graph_legend{detail_index});
        first_test = all_performance.(tests{1});
        first_test_cpu = get_detailed_data(first_test, detail_index);
        title_name = get_result_title(new_test_name);
        figure('Name', title_name, 'Position', [10,10,1280,720]);
        hold on
        for i=1:ntests
            performance = all_performance.(tests{i});
            cpu_performance = get_detailed_data(performance, detail_index);
            plot(cell2mat(performance(:,1)), cpu_performance ./ first_test_cpu, 'color', myColors(i,:));
        end
        xlabel('Input matrix degree');
        ylabel('Speedup (compared to single core)');
        set_figure_props(title_name, tests_names);
    end
end

