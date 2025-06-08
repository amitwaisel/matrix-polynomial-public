function [params] = initialize_params(params)
    if ~isfield(params, 'size')
        params.size = 100;
    end
    if ~isfield(params, 'iterations_to_average')
        params.iterations_to_average = 10;
    end
    if ~isfield(params, 'delta')
        params.delta = 0.1;
    end
    if ~isfield(params, 'nworkers')
        params.nworkers = 28;
    end
    if ~isfield(params, 'degree')
        params.degree = 2000;
    end
    if ~isfield(params, 'coefficients')
        a = -100;
        b = 100;
        params.coefficients = (b-a) .* rand(1,params.degree) + a;
    end
    if ~isfield(params, 'nclusters')
        params.nclusters = params.nworkers;
    end
end

