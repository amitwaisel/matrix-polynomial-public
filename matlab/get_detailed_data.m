function [Y] = get_detailed_data(performance, detail_index)
    n = size(performance, 1);
    X = [];
    Y = [];
        
    for i=1:n
        detailed = cell2mat(performance(i,3));
        if isempty(detailed)
            continue
        end
        index = performance(i,1);
        d = index{1};
        detailed = detailed(detail_index,1);
        detailed = (detailed/1e+9);% / (d^2); % Normalize by n
        X = [X, d];
        Y = [Y, detailed];
    end
    
end