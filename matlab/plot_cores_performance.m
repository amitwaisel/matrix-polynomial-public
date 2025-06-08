function [] = plot_cores_performance(all_performance, test_name)
    tests = fieldnames(all_performance);
    tests_names = cellfun(@get_result_title, tests, 'un', 0);
    ntests = length(tests);
    myColors = jet(ntests);
    first_test = all_performance.(tests{1});
    first_test_cpu = cell2mat(first_test(:,4));
    title_name = get_result_title(test_name);
    figure('Name', title_name, 'Position', [10,10,1280,720]);
    hold on
    for i=1:ntests
        performance = all_performance.(tests{i});
        plot(cell2mat(performance(:,1)), cell2mat(performance(:,4)) ./ first_test_cpu, 'color', myColors(i,:));
    end
    xlabel('Input matrix degree');
    ylabel('Speedup (compared to single core)');
    set_figure_props(title_name, tests_names);
end