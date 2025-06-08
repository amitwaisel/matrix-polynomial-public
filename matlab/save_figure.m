function [] = save_figure(fig_title)
    fig_name = sprintf('figures\\%s - %s.fig', fig_title, datestr(datetime('now')));
    fig_name = strrep(fig_name, ':', '-');
    savefig(fig_name);
    fig_name = sprintf('%s.eps', fig_name);
    print(fig_name,'-depsc2');
end

