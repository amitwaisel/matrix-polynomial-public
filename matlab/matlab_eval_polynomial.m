function [ Q ] = matlab_eval_polynomial( co, A )
%MATLAB_EVAL_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here
    polynomial_string = '0;';
    for i=1:size(co,2)
        polynomial_string = sprintf('%d*A^%d + %s', co(i), i - 1, polynomial_string);
    end
    %fprintf('Evaluating polynomial %s\n', polynomial_string);
    Q = eval(polynomial_string);
end

