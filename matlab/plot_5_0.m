function [] = plot_5_0(params_high, performance_high, matlab_performance_high, params_low, performance_low, matlab_performance_low)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    detailed_plot(performance_high, sprintf('High polynomial degree (%d)', params_high.degree), matlab_performance_high);
    detailed_plot(performance_low, sprintf('Low polynomial degree (%d)', params_low.degree), matlab_performance_low);
end

