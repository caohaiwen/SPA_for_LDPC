% Assumption that our codeword is all-zeros;


% read out LDPC parity-check matrix 
clc;
clear;
%file = fopen("PCMatrix(N=204,K=102,M=102,R=0.5).txt", 'r');
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
        variable(num_row, :) = s;
    end
    if (num_line >=5 && length(s) == col_w)
        num_col = num_col + 1;
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
% %     inv_permutation(num - R_H) = num;
%     inv_permutation(num) = num - R_H;
% end
% 
% while (current_row <= total_row) 
%     if full_H(current_row, current_col) == 1
%         independent_cols(current_row) = current_col ;
% %         inv_permutation(current_row + n - R_H) = current_col;
%         inv_permutation(current_col) = current_row + n -R_H;
%         current_row = current_row +1;
%         current_col = current_col +1;
%     else
%         Q_cols(Q_num) = current_col;
% %         inv_permutation(Q_num) = current_col;
%         inv_permutation(current_col) = Q_num;
%         Q_num = Q_num + 1;
%         current_col = current_col +1;
%     end
% end
% 
% k = n - R_H;  % dimension of our code;
% 
% 
% % construct the permutated generator matrix G_
% G_ = eye(k);
% for num = 1 : k
%     G_(num, k+1:n) = full_H(:,Q_cols(num))';
% end
% 
% % using the inverse permutation to recover the generator matrix
% % corresponding to the original parity-check matrix,i.e. the original code;
% 
% for num = 1 : n
%     G(:,num) = G_(:, inv_permutation(num));
% end    

[G, k] = convertHtoG(H_t);

% AWGN Channel


num = 0;
total_trial = 1000000;
E_b = 1;
R = (n-m)/n ; % designed code rate


for SNR_b = 0: 0.5 : 10   % SNR per information bit,i.e., SNR_b (dB)
    failure = 0;
    num = num + 1;
    for trial = 1 : total_trial
        cm_int = zeros(1,n);
    % randomly generating the information bits;
        u = randi([0,1], 1, k);
        x_s = mod(u*G, 2);
        y_rs = zeros(1, n);
        
    %BPSK modulation
        x_ss = sqrt(R*E_b)*(-2*x_s + 1);    
        
        noise_var = E_b / (2*10^(SNR_b/10));
        for i = 1 : n
            noise = sqrt(noise_var) * randn;
            y_rs(i) = x_ss(i) + noise;
            cm_int(i) = 4*((sqrt(R*E_b))/(2*noise_var))*y_rs(i);
        end
        [decoded_x, convergence] = SPA(cm_int, n, m, row_w, col_w, variable, check);
        if ~(convergence && isequal(decoded_x, x_s))
            failure = failure + 1;
        end
    end
    worderr(num) = failure / total_trial;
end


semilogy(0:0.5:10, worderr, 'rp');

hold on;

title('Word Error Rate vs SNR_b under AWGN Channel using the Sum-Product Algorithm');
xlabel('SNR per information bit (dB)');
ylabel('Word Error Rate');

hold off;

