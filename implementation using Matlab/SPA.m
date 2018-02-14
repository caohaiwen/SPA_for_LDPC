function [decoded_x, convergence] = SPA(cm_int, n, m, row_w, col_w, variable, check)

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
llr_rl = zeros(row, col);
llr_lr = zeros(row, col);
flag = false;

for its = 1 : max_iterations

% massages pass from left(variable nodes) to right(check nodes)
    for k = 1 : row
       variables = variable(k, :); % all the check nodes connected with the k-th variable 
       for j = 1 : degree_variable % compute the message from variable node k to its j_th check node
           llr_lr(k, variables(j)) = cm_int(k) + llr_rl(k, variables(mod(j+1-1, degree_variable)+1));
           for t = j+2 : j+degree_variable-1
               llr_lr(k, variables(j)) = llr_lr(k, variables(j)) + llr_rl(k, variables(mod(t-1, degree_variable)+1));
           end
       end 
    end
    
% massages pass from right(check nodes to left(variable nodes)
    for j = 1 : col
        checks = check(j, :);
        for k = 1 : degree_check % compute the message from check node j to its k_th variable node
            llr_rl(checks(k), j) = tanh(llr_lr(checks(mod(k+1-1, degree_check)+1), j)/2);
            for t = k+2 : k+degree_check-1
                llr_rl(checks(k), j) = llr_rl(checks(k), j) * tanh(llr_lr(checks(mod(t-1, degree_check)+1), j)/2);
            end
            llr_rl(checks(k), j ) = 2*atanh(llr_rl(checks(k), j));
        end     
    end
    
% compute the current beliefs in each varaible node 
    for i = 1 : row
        variables = variable(i, :);
        beliefs(i) = cm_int(i);
        for m = 1 : degree_variable % compute the i_th variable node's belief 
            beliefs(i) = beliefs(i) + llr_rl(i, variables(m));
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