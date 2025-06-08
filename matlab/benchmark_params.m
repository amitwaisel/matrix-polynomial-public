function [performance] = benchmark_params(blocks_algorithm, sylvester, schur, eigenvalues_sort, delta, coefficients, nworkers, params )
    params.blocks_algorithm = blocks_algorithm;
    params.sylvester = sylvester;
    params.schur = schur;
    params.eigenvalues_sort = eigenvalues_sort;
    params.delta = delta;
    params.coefficients = coefficients;
    params.nworkers = nworkers;
    performance = benchmark(params);
end