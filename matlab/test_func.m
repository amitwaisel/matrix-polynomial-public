function [] = test_func(test_name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if exist('test_name', 'var') && ~isempty(test_name)
        fprintf('test_name is %s', test_name);
    else
        fprintf('No test_name given');
    end
end

