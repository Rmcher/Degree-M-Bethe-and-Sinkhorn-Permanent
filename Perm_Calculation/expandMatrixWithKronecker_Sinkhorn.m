function expandedMatrix = expandMatrixWithKronecker_Sinkhorn(A, n, M)
% expandMatrixWithKronecker_Sinkhorn - Expands matrix A to dimension nM *
% nM for M cover Sinkhorn approximation
    %
    % Syntax: expandedMatrix = expandMatrixWithKronecker_Sinkhorn(A, n, M)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %   M - The expand dimension, cover of graphs.
    %
    % Output:
    %   expandedMatrix - expanded matrix using Kronecker product 1/M allone
    %   matrix.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.15
    
    % Create an all-one matrix with elements equal to 1/M
    allOneMatrix = ones(M) / M;

    % Initialize the expanded matrix
    expandedMatrix = zeros(n * M);

    % Iterate through each element of A to create the expanded matrix
    for i = 1:n
        for j = 1:n
            % Perform Kronecker product between A(i, j) and the all-one matrix
            kronProduct = A(i, j) * allOneMatrix;

            % Place Kronecker product outcome to the new matrix
            rowStart = (i - 1) * M + 1;
            rowEnd = i * M;
            colStart = (j - 1) * M + 1;
            colEnd = j * M;

            expandedMatrix(rowStart:rowEnd, colStart:colEnd) = kronProduct;
        end
    end

    return;
end