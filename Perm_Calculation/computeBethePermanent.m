function [permB, V]  = computeBethePermanent(A, n)
% computeBethePermanent - Computes the Bethe approximation permanent of
% matrix A based on Sum Product Algorithm (SPA) in Factor Graphs.
    %
    % Syntax: [permB, V] = computeBethePermanent(A, n)
    %
    % Inputs:
    %   A - The input matrix (n x n).
    %   n - The dimension of the matrix A (integer).
    %
    % Output:
    %   permB - The computed permanent using partition function in Factor
    %   Graphs.
    %
    % Author: WU Binghong
    % Date: 2024.Oct.02
    
    %--------------------------------------------------------------------
    % Initialize the n*n message cell and V symbol matrix, cell{i, j} 
    % contain the message on x(i, j), in the form of vector = [u(0), u(1)]T
    MessageCellRight = cell(n, n);
    MessageCellLeft = cell(n, n);
    VRight = {};
    VLeft = {};
    
    % Store values for plotting
    V = [];

    % Random Initialize the message on each edge for t=0, sending right,
    % i.e. x, ensuring sum=1. For u_{ij}(x_ij) = 1, MessageCell{i, j}(2)
    for i = 1:n
        for j = 1:n
            % r = rand(); 
            r = 1/3; 
            vector = [r; 1 - r]; 
            MessageCellRight{i, j} = vector;
            VRight{1,1}(i, j) = r/(1-r);
        end
    end 

    %--------------------------------------------------------------------
    % Message sending algorithm in iteration
    
    % Iteration setup
    t = 0;  % Initial time step
    tolerance = 1e-4;  % Stopping condition tolerance for VLeft comparison

    while true
        t = t + 1;

        if mod(t, 2) == 1  % t is odd: update MessageCellLeft and VLeft
            % Update the MessageCellLeft and VLeft based on current MessageCellRight
            for i = 1:n
                for j = 1:n
                    % Initialize MessageCellLeft{i,j}(1)
                    sumTerm = 0;
                    for i_prime = 1:n
                        if i_prime ~= i && A(i_prime, j) ~= 0
                            % Calculate \prod(MessageCellRight{i'',j}(1)), i'' != i or i'
                            productTerm = 1;
                            for i_double_prime = 1:n
                                if i_double_prime ~= i && i_double_prime ~= i_prime && A(i_double_prime, j) ~= 0
                                    productTerm = productTerm * MessageCellRight{i_double_prime, j}(1);
                                end
                            end
                            % ..times MessageCellRight{i',j}(2) and sqrt(A(i',j))
                            productTerm = productTerm * MessageCellRight{i_prime, j}(2) * sqrt(A(i_prime, j));
                            % Iterate sumTerm
                            sumTerm = sumTerm + productTerm;
                        end
                    end
                    % Update MessageCellLeft{i,j}(1)
                    MessageCellLeft{i,j}(1) = sumTerm;
        
                    % Calculate MessageCellLeft{i,j}(2)
                    % ..times sqrt(A(i,j))
                    productTerm = sqrt(A(i,j));
                    % Calculate \prod(MessageCellRight{i',j}(1)), i' != i 
                    for i_prime = 1:n
                        if i_prime ~= i && A(i_prime, j) ~= 0
                            productTerm = productTerm * MessageCellRight{i_prime, j}(1);
                        end
                    end
                    MessageCellLeft{i,j}(2) = productTerm;
        
                    % Normalize MessageCellLeft{i,j}
                    total = MessageCellLeft{i,j}(1) + MessageCellLeft{i,j}(2);
                    MessageCellLeft{i,j} = MessageCellLeft{i,j} / total;
        
                    % Update VLeft
                    VLeft{(t + 1) / 2,1}(i,j) = MessageCellLeft{i,j}(1) / MessageCellLeft{i,j}(2);
                end
            end

            if nargout > 1
                % Calculate V = VLeft * VRight
                productMatrix = VLeft{(t + 1) / 2} .* VRight{(t + 1) / 2};
                % Store V for plot
                V = [V; productMatrix(:)'];
            end

        else  % t is even: update MessageCellRight and VRight
            % Update the MessageCellRight and VRight based on current MessageCellLeft
            for j = 1:n
                for i = 1:n
                    % Initialize MessageCellRight{i,j}(1)
                    sumTerm = 0;
                    for j_prime = 1:n
                        if j_prime ~= j && A(i, j_prime) ~= 0
                            % Calculate \prod(MessageCellRight{i,j''}(1)), j'' != j or j'
                            productTerm = 1;
                            for j_double_prime = 1:n
                                if j_double_prime ~= j && j_double_prime ~= j_prime && A(i, j_double_prime) ~= 0
                                    productTerm = productTerm * MessageCellLeft{i, j_double_prime}(1);
                                end
                            end
                            % ..times MessageCellLeft{i,j'}(2) and sqrt(A(i,j'))
                            productTerm = productTerm * MessageCellLeft{i, j_prime}(2) * sqrt(A(i, j_prime));
                            % Iterate sumTerm
                            sumTerm = sumTerm + productTerm;
                        end
                    end
                    % Update MessageCellRight{i,j}(1)
                    MessageCellRight{i,j}(1) = sumTerm;
        
                    % Calculate MessageCellRight{i,j}(2)
                    % ..times sqrt(A(i,j))
                    productTerm = sqrt(A(i,j));
                    % Calculate \prod(MessageCellRight{i,j'}(1)), j' != i 
                    for j_prime = 1:n
                        if j_prime ~= j && A(i, j_prime) ~= 0
                            productTerm = productTerm * MessageCellLeft{i, j_prime}(1);
                        end
                    end
                    MessageCellRight{i,j}(2) = productTerm;
        
                    % Normalize MessageCellRight{i,j}
                    total = MessageCellRight{i,j}(1) + MessageCellRight{i,j}(2);
                    MessageCellRight{i,j} = MessageCellRight{i,j} / total;
        
                    % Update VRight
                    VRight{t / 2 + 1,1}(i,j) = MessageCellRight{i,j}(1) / MessageCellRight{i,j}(2);
                end
            end
        end

        %----------------------------------------------------------------
        % Check stopping condition
        if t > 2 && mod(t, 2) == 1  % Only check VLeft in odd iterations
            diff = max(max(abs(VLeft{(t + 1) / 2} - VLeft{(t - 1) / 2})));  % Difference between VLeft at t and t-2
            if diff < tolerance || t > 300000
                % disp(['Iteration converged at t = ', num2str(t)]);
                break;  % Stop if the change is smaller than tolerance
            end
        end
    end

    %--------------------------------------------------------------------
    % Calculate permB = Zf/Ze
    Zf_left = 1;  % Initialize Zf_left
    Zf_right = 1; % Initialize Zf_right
    % Zf_left
    for i = 1:n
        sum_i = 0;
        for j = 1:n
            term = sqrt(A(i,j)) * MessageCellLeft{i,j}(2);
            product_j_prime = 1;
            for j_prime = 1:n
                if j_prime ~= j && A(i, j_prime) ~= 0
                    product_j_prime = product_j_prime * MessageCellLeft{i,j_prime}(1);
                end
            end
            sum_i = sum_i + term * product_j_prime;
        end
        Zf_left = Zf_left * sum_i;
    end

    % Zf_right
    for j = 1:n
        sum_j = 0;
        for i = 1:n
            term = sqrt(A(i,j)) * MessageCellRight{i,j}(2);
            product_i_prime = 1;
            for i_prime = 1:n
                if i_prime ~= i && A(i_prime, j) ~= 0
                    product_i_prime = product_i_prime * MessageCellRight{i_prime,j}(1);
                end
            end
            sum_j = sum_j + term * product_i_prime;
        end
        Zf_right = Zf_right * sum_j;
    end

    Ze = 1;
    for i = 1 : n
        for j = 1 : n
            Ze = Ze * (MessageCellRight{i,j}(1)* MessageCellLeft{i,j}(1)...
                + MessageCellRight{i,j}(2)* MessageCellLeft{i,j}(2));
        end
    end
    
    Zf = Zf_left * Zf_right;

    permB = Zf/Ze;

    return;
end
