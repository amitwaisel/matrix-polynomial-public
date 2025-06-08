function [ Q ] = matlab_polynomial( co, A )
%MATLAB_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here
    Q = zeros(size(A));
    temp = eye(size(A));
    for i=1:size(co,2)
        Q = Q + temp * co(i);
        temp = temp * A;
    end
end

