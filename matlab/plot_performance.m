function [] = plot_performance(all_performance, test_name, matlab_performance, options)
    tests = fieldnames(all_performance);
    tests_names = cellfun(@get_result_title, tests, 'un', 0);
    ntests = length(tests);
    myColors = jet(ntests);
    
    if ~exist('options','var')
        options = struct();
    end
    
    title_name = get_result_title(test_name);
    figure('Name', title_name, 'Position', [10,10,1280,720]);
    hold on
    for i=1:ntests
        performance = all_performance.(tests{i});
        plot(cell2mat(performance(:,1)), cell2mat(performance(:,4)), 'color', myColors(i,:));
    end
    xlabel('dimension');
    ylabel('time (s)');
    set_figure_props(title_name, tests_names);
    %figure
    for i=1:ntests
        %subplot(ntests/2, ntests-ntests/2, i)
        current_name = tests{i};
%         if exist('test_name', 'var') && ~isempty(test_name)
%             current_name = test_name;
%         end
        if exist('matlab_performance', 'var')
            detailed_plot(all_performance.(tests{i}), current_name, matlab_performance, options);
        else
            detailed_plot(all_performance.(tests{i}), current_name, [], options);
        end
    end
        
end