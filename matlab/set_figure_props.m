function [] = set_figure_props(title_name, graph_legend, do_linear)
    set(gca, 'FontSize', 20); 
    if exist('do_linear','var') && do_linear
        set(gca, 'Yscale', 'linear');
    else
        set(gca, 'Yscale', 'log');
    end        
    h=legend(graph_legend, 'Location', 'southeast');
    set(h,'FontSize',16);
    title(title_name, 'Interpreter', 'none');
    save_figure(title_name);
end

