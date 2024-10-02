clear;
clc;

% Define a specified square matrix
A = [1, 1, 1; 1, 1, 1; 1, 1, 1];
% Get the matrix size
[n, ~] = size(A);

% Calculate the perm in normal form and Bethe Approximation
perm = computeNormalPermanent(A);
[permB, V] = computeBethePermanent(A, n);

% Plot the VLeft values over time for each element
plotVLeft(V, n);

% Print the outcome
fprintf('Perm(A) = %.2f\n', perm);
fprintf('PermB(A) = %.2f\n', permB);