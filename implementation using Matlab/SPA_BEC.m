function [decoded_x, convergence] = SPA_BEC(y_r,n, m, row_w, col_w, variable, check, max_iterations, P)

% decoded_x: the output of the decoder
% convergence : whether it can convergent to a codeword within the max_iterations;
% cm_int : the intrisic information from the channel model.
% y_r : the received signal
% 
%SPA implementation 

row = n;
col = m;
degree_variable = row_w;
degree_check = col_w;

decoded_x = ones(1, row);
convergence = false;
llr_lr = zeros(1, row*degree_variable);
llr_rl = zeros(1, col*degree_check);

edge = 0;
next_edge = 0;
cnt = 0; 
prv_cnt = 0;

codeword = y_r;

% iteration0 : from left to right is just the message from the channel. 
for k = 1: row 
    variables = variable(k, :);
    for j = 1 : degree_variable
        if (variables(j) ~= 0)
          edge = P(k, variables(j));
          llr_lr(edge) = codeword(k);
        end
    end
end

for its = 1 : max_iterations
    
% massages pass from right(check nodes to left(variable nodes)
    for j = 1 : col
        checks = check(j, :);
        for k = 1 : degree_check % compute the message from check node j to its k_th variable node
            if (checks(k) ~= 0)
            % the outgoing message from the check node is an erasure if there
            % exists any erasure from the incoming message, otherwise, all the 
            % incoming messages are not erased and then mod(all the sum of 
            % incoming messages).
                edge = P(checks(k), j);
                llr_rl(edge) = 0;
                
                for t = k+1 : k+degree_check-1
                    
                    if (checks(mod(t-1, degree_check)+1) == 0)
                        continue;
                    end
                    next_edge = P(checks(mod(t-1, degree_check)+1), j);

                    if (llr_lr(next_edge) ~= 2)
                        llr_rl(edge) = llr_rl(edge) + llr_lr(next_edge);
                        llr_rl(edge) = mod(llr_rl(edge), 2);
                    else 
                        llr_rl(edge) = 2;
                        break;
                    end
                    
                end
            end     
        end
    end
    
% massages pass from left(variable nodes) to right(check nodes)
    for k = 1 : row
      variables = variable(k, :);
       
      % the outgoing message from the variable node is an erasure if all
      % incoming messages are erasures, thus we initiate every outgoing
      % message an erasure at first.
      for j = 1 : degree_variable % compute the message from variable node k to its j_th check node
          if (variables(j) ~= 0)
            edge = P(k, variables(j));
            if (codeword(k) ~= 2)
              llr_lr(edge) = codeword(k);
            else 
              for t = j+1 : j+degree_variable-1
                if (variables(mod(t-1, degree_variable)+1) == 0)
                  continue;
                end
                next_edge = P(k, variables(mod(t-1, degree_variable)+1));

                if (llr_rl(next_edge) ~= 2)
                  llr_lr(edge) = llr_rl(next_edge);
                  codeword(k) =  llr_lr(edge);
                  break;
                end
              end
            end
          end

      end
      if (codeword(k) ~= 2)
          cnt = cnt + 1;
      end
    end
    
      
% check whether this codeword is a real codeword or not
    if (cnt == row)    
        decoded_x = codeword;
        convergence = true;
        break;
    else if (cnt == prv_cnt)
            break;
        end
    end
    
    prv_cnt = cnt;
    cnt = 0;
end

end