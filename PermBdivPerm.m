% author: WU Binghong
% Date: 2024.Oct.02

% Define the number of matrices
numMatrices = 1000;
n = 5;

% Define the arrays to store resultsmatrixSize
permValues = zeros(1, numMatrices);
permBValues = zeros(1, numMatrices);

% Generate 1000 random matrices and calculate relevant results
for i = 1:numMatrices
    % Generate a 5x5 matrix with entries from uniform distribution [0, 1]
    A = rand(n, n);
    
    % Compute the permanent of A
    permValues(i) = computeNormalPermanent(A);
    
    % Compute Bethe permanent approximations of A
    permBValues(i) = computeBethePermanent(A, n);
end

% Plot the results
figure;
hold on;

% Plot perm(A) vs permB(A) (cyan triangles)
loglog(permValues, permBValues, 'c^', 'DisplayName', 'permB(A)');

% Theoretical ratio for perm/permB
gammaB = sqrt(2 * pi * n) / exp(1);

% Plot theoretical lines
loglog([min(permValues), max(permValues)], ...
    [min(permValues)/gammaB, max(permValues)/gammaB], 'k--', 'DisplayName', 'Theoretical perm/permB');

% Set labels and legend
xlim([10^-(1.3) 10^(1.3)]);
ylim([10^-(1.3) 10^(1.3)]);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('perm(A)');
ylabel('permB(A) and !permB,2(A)');
legend;
title('Comparison of perm(A) vs permB(A) and !permB,2(A)');

hold off;


