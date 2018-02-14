function [G, k] = convertHtoG(H)
% this function is to convert a general parity-check matrix H into a
% corresponding generator matrix G
% input : a general parity-check matrix H
% output : the corresponding generator matrix G, the dimension k

H_ = H_sf(H);  

current_row = 0;
[num_row, num_col] = size(H);
n = num_col;

% removing all the all-zero rows in H_
for row = 1 : num_row
    if (sum(H_(row,:)) ~= 0)
        current_row = current_row + 1;
        full_H(current_row,:) = H_(row,:);
    end
end

total_row = current_row; 

% the rank of parity-check matrix H
R_H = total_row;

% select all the independent columns from full-row-rank full_H
% current_col = 1;
% current_row = 1;
% Q_num = 1;
independent_cols =zeros(1,total_row);
% Q_cols = zeros(1, n-R_H);
inv_permutation = zeros(1, n); 

% for num = R_H + 1 : n
%     Q_cols(num-R_H) = num;
% %     inv_permutation(num - R_H) = num;
%     inv_permutation(num) = num - R_H;
% end

Q_cols = 1 : n;
for num = 1 : R_H
    if full_H(num, num) == 1
        independent_cols(num) = num;
        inv_permutation(num) = num + n - R_H;
        location = find(Q_cols == num);
        Q_cols(location) = [];   %delete this independent column from Q_cols
    else 
        location = 0;
        for col = Q_cols
            if full_H(num, col) == 1
                independent_cols(num) = col;
                inv_permutation(col) = num + n -R_H;
                location = find(Q_cols == col);
                break;
            end
        end
        if location ~= 0
            Q_cols(location) = [];
        end
    end
end

col = 0;
for num = Q_cols
    col = col +1;
    inv_permutation(num) = col;
end
        
% while (current_row <= total_row) 
%     
%     if full_H(current_row, current_col) == 1
%         independent_cols(current_row) = current_col ;
%         inv_permutation(current_col) = current_row + n -R_H;
%         current_row = current_row +1;
%         current_col = current_col +1;
%     else
%         Q_cols(Q_num) = current_col;
%         inv_permutation(current_col) = Q_num;
%         Q_num = Q_num + 1;
%         current_col = current_col +1;
%     end
% end

k = n - R_H;  % dimension of our code;


% construct the permutated generator matrix G_
G_ = eye(k);
for num = 1 : k
    G_(num, k+1:n) = full_H(:,Q_cols(num))';
end

% using the inverse permutation to recover the generator matrix
% corresponding to the original parity-check matrix,i.e. the original code;

for num = 1 : n
    G(:,num) = G_(:, inv_permutation(num));
end    
