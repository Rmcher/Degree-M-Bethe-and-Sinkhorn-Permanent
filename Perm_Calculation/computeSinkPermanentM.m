function permSM = computeSinkPermanentM(A, n, M)
% computeSinkPermanentM - Computes the Sinkhorn approximation permanent of
% matrix A.
    %
    % Syntax: permSM = computeSinkPermanentM(A, n, M)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %   M - Expand dimension, M cover of graphs.
    %
    % Output:
    %   permSM - The M cover Sinkhorn Approximation of Permenent
    %
    % Author: WU Binghong
    % Date: 2024.Oct.15

    expandedMatrix = expandMatrixWithKronecker_Sinkhorn(A, n, M);
    
    % Compute perm of expandedMatrix in standard way
    expandedSize = n * M;
    standardperm = perm(expandedMatrix, expandedSize);
    
    % Calculate permSM
    permSM = standardperm ^ (1 / M);

    return;
end