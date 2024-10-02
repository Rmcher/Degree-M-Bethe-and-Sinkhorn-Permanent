% Author: WU Binghong
% Date: 2024.Oct.02

% Given the V information and the matrix size, plot the convergence trend of V
function plotV(V, n)
    figure;
    hold on;
    
    % Get the time of iteration
    numIterations = size(V, 1);

    for idx = 1:(n * n)
        plot(1:numIterations, V(:, idx), 'LineWidth', 2);
    end

    xlabel('Iteration (t)');
    ylabel('V');
    title('V Over Time');

    % Unicode
    legendLabels = {};
    for i = 1:n
        for j = 1:n
            legendLabels{end+1} = ['\bf\mu_{', num2str(i), ',', num2str(j), '}'];
        end
    end

    legend(legendLabels, 'Location', 'bestoutside', 'Interpreter', 'tex');

    hold off;
end
