function [decoded_x, convergence] = SPA(cm_int, n, m, row_w, col_w, variable, check, P)

% decoded_x: the output of the decoder
% convergence : whether it can convergent to a codeword within the max_iterations;
% cm_int : the intrisic information from the channel model.
% y_r : the received signal
% 
%SPA implementation 

max_iterations = 100;
row = n;
col = m;
degree_variable = row_w;
degree_check = col_w;

decoded_x = ones(1, row);
convergence = false;

beliefs = zeros(1, row);
llr_rl = zeros(1, col*degree_check);
llr_lr = zeros(1, row*degree_variable);
flag = false;
open_variableN = zeros(1, 2*(degree_variable - 1));
open_checkN = zeros(1, 2*(degree_check - 1));

for its = 1 : max_iterations

% massages pass from left(variable nodes) to right(check nodes) using "opening the box" operation 
    for k = 1 : row
        variables = variable(k, :); % all the check nodes connected with the k-th variable 

        % compute all the messages (after opening the box) from top to bottom
        for j = 1 : degree_variable - 1 
            if (variables(j) == 0) 
                continue;
            end

            edge = P(k, variables(j));
            if (j == 1)
                open_variableN(j) = cm_int(k) + llr_rl(edge);
            else
                open_variableN(j) = open_variableN(j-1) + llr_rl(edge);
            end
        end

        % compute all the messages (after opening the box) from bottom upto top 
        num = degree_variable ;
        for j = degree_variable : -1 : 2
            if (variables(j) == 0) 
                continue;
            end
            edge = P(k, variables(j));
            % messages from bottom upto top, and merge the upper and lower messages to compute the outgoing messages;
            if (j == degree_variable)
                open_variableN(num) = llr_rl(edge);
                llr_lr(edge) = open_variableN(degree_variable - 1);
            else
                open_variableN(num) = open_variableN(num - 1) + llr_rl(edge);
                llr_lr(edge) = open_variableN(num - 1) + open_variableN(2*(degree_variable -1) - (num - 1));
            end
            num = num + 1;
        end

        % compute the first outgoing message
        if (variables(1) == 0)
            continue;
        end
        edge = P(k, variables(1));
        llr_lr(edge) = open_variableN(num - 1) + cm_int(k);
    end
    

% massages pass from right(check nodes to left(variable nodes) using "opening the box" operation
    for j = 1 : col
        checks = check(j, :);

        % open the box : messages passing from top to bottom
        for k = 1 : degree_check - 1 
            if (checks(k) == 0)
                continue;
            end
            edge = P(checks(k), j);
            value = llr_lr(edge) / 2.0;
            tanhVaule = 1.0;
            if (abs(value) < 17.5)  % otherwise, tanh(tanhVaule) = +1/-1
                tanhVaule = tanh(value);
            elseif (abs(value) < -17.5)
                tanhVaule = -1.0;
            end

            if (k == 1)
                open_checkN(k) = tanhVaule ;
            else
                open_checkN(k) = open_checkN(k -1) * tanhVaule;
            end
        end

        %open the box: messages passing from bottom to top, and also compute all the outgoing messages;
        num = degree_check; 
        for k = degree_check : -1 : 2
            if (checks(k) == 0)
                continue;
            end
            edge = P(checks(k), j);
            value = llr_lr(edge) / 2.0;
            tanhVaule = 1.0;
            if (abs(value) < 17.5)  % otherwise, tanh(tanhVaule) = +1/-1
                tanhVaule = tanh(value);
            elseif (abs(value) < -17.5)
                tanhVaule = -1.0;
            end

            if (k == degree_check)
                open_checkN(num) = tanhVaule;
                tanhVaule = open_checkN(degree_check - 1);
            else
                open_checkN(num) = open_checkN(num -1) * tanhVaule;
                tanhVaule = open_checkN(num - 1) * open_checkN(2*(degree_check -1) - (num - 1));
            end
            
            if (abs(tanhVaule - 1.0) <= 1e-15)
                value = 17.5;
            elseif (abs(tanhVaule - (-1.0)) <= 1e-15)
                value = -17.5;
            else
                value = atanh(tanhVaule);
            end
            llr_rl(edge) = 2 * value;
            num = num + 1;
        end

        % compute the first outgoing message
        if (checks(1) == 0)
            continue;
        end
        edge = P(checks(1), j);
        tanhVaule = open_checkN(num - 1);
        if (abs(tanhVaule - 1.0) <= 1e-15)
            value = 17.5;
        elseif (abs(tanhVaule - (-1.0)) <= 1e-15)
            value = -17.5;
        else
            value = atanh(tanhVaule);
        end
        llr_rl(edge) = 2 * value;
    end
    
% compute the current beliefs in each varaible node 
    for i = 1 : row
        variables = variable(i, :);
        beliefs(i) = cm_int(i);
        for m = 1 : degree_variable  % compute the i_th variable node's belief 
            if (variables(m) == 0)
                continue;
            end
            edge = P(i, variables(m));
            beliefs(i) = beliefs(i) + llr_rl(edge);
        end
    end
    
    codeword = ones(1, row);
    for i = 1 : row
        if beliefs(i) >= 0
            codeword(i) = 0;
        else
            codeword(i) = 1;
        end
    end
      
% check whether this word is a real codeword or not
    for j = 1 : col
        checks = check(j, :);
        sum_checks = 0;
        for k = 1 : degree_check
            sum_checks = sum_checks + codeword(checks(k));
        end
        if (mod(sum_checks, 2)==1)
            flag = true;
            break;
        end
    end
    
    if (flag == false)
        decoded_x = codeword;
        convergence = true;
        break;
    end
    
    flag = false;
end

end