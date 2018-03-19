% Assumption that our codeword is all-zeros;


% read out LDPC parity-check matrix 
clc;
clear;
file = fopen("PCMatrix(N=204,K=102,M=102,R=0.5).txt", 'r');
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
P = zeros(n, m);
edge = 0;

for i = 1 : n
    for j = 1 : m
        if ismember(j, variable(i,:))
            H(i, j) = 1;
            edge = edge + 1;
            P(i, j) = edge;
        end
    end
end

H_t = H';

[G, k] = convertHtoG(H_t);


% BSC: cross prob. = 0.05: 0.05 : 0.40;
u = randi([0,1], 1, k);

x_s = mod(u*G, 2);


cm_int = zeros(1,n);
total_trial = 1000;
num = 0;
for cross_prob = 0.05: 0.02 : 0.25
    failure = 0;
    num = num + 1;
    for trial = 1 : total_trial
        y_r = bsc(x_s, cross_prob);
        for i = 1 : n
            if (y_r(i) == 1)
                cm_int(i) = log(cross_prob / (1 - cross_prob));
            else
                cm_int(i) = log((1 - cross_prob) / cross_prob);
            end
        end
        [decoded_x, convergence] = SPA(cm_int, n, m, row_w, col_w, variable, check, P);
        if ~(convergence && isequal(decoded_x, x_s))
            failure = failure + 1;
        end
    end
    worderr(num) = failure/total_trial;
end


plot(0.05: 0.02 : 0.25, worderr);


