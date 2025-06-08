function [ B ] = build_unit_eigenvalue_matrix( n, nclusters )
    % For each block [a,b;c,a], eigenvalues are a+-i*sqrt(b*c).
    % radius=1, b*c<0
    % norm(a+i*sqrt(bc))=1 -> a^2+2ai*sqrt(bc)-bc=1 ->
    % bc-2ai*sqrt(bc)=a^2-1 -> 
    % a = cos(theta), b = rand, c = 1/[b*(sin(theta)^2)].
    if n < nclusters
        error('Cannot generate %d clusters for %d elements', nclusters, n);
    end
    angles = 2*pi*(0:1/nclusters:1) + 0.5/nclusters;
    angles(end) = [];
    points_complex = exp(i*angles);
    points_a = real(points_complex); % cos(angles); % ones(size(angles))./cos(angles);
    points_imag = imag(points_complex); % c = imag^2/b
    points_b = rand(size(points_a));
    % points_c = ones(size(angles)) ./ (points_b .* (sin(angles)).^2);
    points_c = (points_imag.^2) ./ points_b;
    
    nblocks = floor(n/2);
    blocks_in_each_cluster = floor(nblocks/nclusters);
    blocks_in_each_cluster_vector = ones(nclusters,1)*blocks_in_each_cluster;
    blocks_in_each_cluster_vector(end) = nblocks - blocks_in_each_cluster*(nclusters-1);
    
    a = repelem(points_a, blocks_in_each_cluster_vector);
    b = repelem(points_b, blocks_in_each_cluster_vector);
    c = repelem(points_c, blocks_in_each_cluster_vector);
    
    perm = randperm(nblocks);
    a=a(perm); b=-1*b(perm); c=c(perm);
    
    A = rand(n);
    [Q,T] = schur(A, 'real');
    
    T = 1e-5 * triu(ones(n, n));
    for j=1:2:n
        k=(j+1)/2;
        T(j:j+1,j:j+1) = [a(k), b(k); c(k), a(k)];
    end
    
    if nblocks*2<n % last element not in 2x2 block
        T(end,end)=1;
    end
    %T = T / norm(T,2); %normalize_triangular_matrix(T);
    %T = normalize_triangular_matrix(T);
    
    B = Q*T*Q';
    return;
    
    close all;
    figure
    plot(eig(T),'x')
    hold on;
    plot(eig(B),'o')
    return;
    
    d = eye(n) .* T;
    d = d + (eye(n) == 0);
    T = T ./ sqrt(d .* d);
end

function [ nT ] = normalize_triangular_matrix( T )
%NORMALIZE_TRIANGULAR_MATRIX Normalize all blocks on the diagonal
%   Calculate the norm of every 1x1/2x2 block on the diagonal and normalize
    in_block = false;
    nT = T;
    for i=1:size(T,1)
        if in_block
            in_block = false;
            continue;
        end
        block_start = i;
        block_end = i;
        if i < size(T,1) && T(i+1,i) ~= 0
            block_end = i + 1;
            in_block = true;
        end
        nT(block_start:block_end,block_start:block_end) = T(block_start:block_end,block_start:block_end) / norm(T(block_start:block_end,block_start:block_end));
    end
end