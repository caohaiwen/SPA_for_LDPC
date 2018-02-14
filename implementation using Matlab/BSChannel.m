% Assumption that our codeword is all-zeros;


% read out LDPC parity-check matrix 
clc;
clear;
file = fopen("PCMatrix(96.3.963 (N=96,K=48,M=48,R=0.5)).txt", 'r');
num_line = 0;
num_row = 0;
num_col = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
    if (num_line == 1)
        n = s(1);
        m = s(2);
    end
    if (num_line == 2)
        row_w = s(1);
        col_w = s(2);
        variable = zeros(n,row_w);
        check = zeros(m, col_w);
    end
    if (num_line >= 5 &&  length(s) == row_w)
        num_row = num_row + 1;
    % the (num_row)-th variable node with its connected (row_w) check nodes
        variable(num_row, :) = s;  
    end
    if (num_line >=5 && length(s) == col_w) 
        num_col = num_col + 1;
    % the (num_col)-th check node with its connected (col_w) variable nodes 
        check(num_col, :) = s;
    end           
end

fclose(file);

% construct the parity-check matrix 

H = zeros(n, m);

for i = 1 : n
    for j = 1 : m
        if ismember(j, variable(i,:))
            H(i, j) = 1;
        end
    end
end

H_t = H';

[G, k] = convertHtoG(H_t);

% H_ = H_sf(H_t);  
% 
% current_row = 0;
% 
% % removing all the all-zero rows in H_
% for row = 1 : num_col
%     if (sum(H_(row,:)) ~= 0)
%         current_row = current_row + 1;
%         full_H(current_row,:) = H_(row,:);
%     end
% end
% 
% total_row = current_row; 
% 
% % the rank of parity-check matrix H
% R_H = total_row;
% 
% % select all the independent columns from full-row-rank full_H
% current_col = 1;
% current_row = 1;
% Q_num = 1;
% independent_cols =zeros(1,total_row);
% Q_cols = zeros(1, n-R_H);
% inv_permutation = zeros(1, n); 
% 
% for num = R_H + 1 : n
%     Q_cols(num-R_H) = num;
%     %permutation(num - R_H) = num;
%     inv_permutation(num) = num - R_H;
% end
% 
% while (current_row <= total_row) 
%     if full_H(current_row, current_col) == 1
%         independent_cols(current_row) = current_col ;
%       %  permutation(current_row + n - R_H) = current_col;
%         inv_permutation(current_col) = current_row + n -R_H;
%         current_row = current_row +1;
%         current_col = current_col +1;
%     else
%         Q_cols(Q_num) = current_col;
%       %  permutation(Q_num) = current_col;
%         inv_permutation(current_col) = Q_num;
%         Q_num = Q_num + 1;
%         current_col = current_col +1;
%     end
% end
% 
% 
% 
% k = n - R_H;  % dimension of our code;
% 
% 
% % construct the permutated generator matrix G_p
% G_p = eye(k);
% for num = 1 : k
%     G_p(num, k+1:n) = full_H(:,Q_cols(num))';
% end
% 
% % using the inverse permutation to recover the generator matrix
% % corresponding to the original parity-check matrix,i.e. the original code;
% 
% for num = 1 : n
%     G(:,num) = G_p(:, inv_permutation(num));
% end 


% BSC: cross prob. = 0.05: 0.05 : 0.40;

total_trial = 1000000;
num = 0;

for cross_prob = 0.05: 0.02 : 0.45
    failure = 0;
    num = num + 1;
    for trial = 1 : total_trial
        cm_int = zeros(1,n);
        
        u = randi([0,1], 1, k);
        x_s = mod(u*G, 2);
        y_r = bsc(x_s, cross_prob);
        for i = 1 : n
            if (y_r(i) == 1)
                cm_int(i) = log(cross_prob / (1 - cross_prob));
            else
                cm_int(i) = log((1 - cross_prob) / cross_prob);
            end
        end
        [decoded_x, convergence] = SPA(cm_int, n, m, row_w, col_w, variable, check);
        if ~(convergence && isequal(decoded_x, x_s))
            failure = failure + 1;
        end
    end
    worderr(num) = failure/total_trial;
end



semilogy(0.05: 0.02 : 0.45, worderr, 'rp');

hold on; 

title('Word Error Rate vs crossover probabilty under BSC Channel using the Sum-Product Algorithm');
xlabel('SNR per information bit (dB)');
ylabel('Word Error Rate');

hold off;



