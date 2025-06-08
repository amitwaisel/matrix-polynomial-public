function [] = plot_ps_performance(all_performance)
    tests = fieldnames(all_performance);
    tests_names = cellfun(@get_result_title, tests, 'un', 0);
    ntests = length(tests);
    myColors = jet(ntests);
    title_name = 'Paterson-Stockmeyer implementations';
    figure('Name', title_name, 'Position', [10,10,1280,720]);
    xlabel('dimension');
    ylabel('time (ms)');
    hold on
    ps_inedx = 6;
    for i=1:ntests
        performance = all_performance.(tests{i});
        stats = cell2mat(performance(:,2)); % get_detailed_data(performance, ps_inedx);
        plot(cell2mat(performance(:,1)), stats, 'color', myColors(i,:), 'LineWidth', 4*ntests-4*(i-1));
    end
    set_figure_props(title_name, tests_names);
end
