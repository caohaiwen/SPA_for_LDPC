% Assumption that our codeword is all-zeros;


% read out LDPC parity-check matrix 
clc;
clear;
file = fopen(".\Data (i.e. LDPC Matrices)\PCMatrix(96.3.963 (N=96,K=48,M=48,R=0.5)).txt", 'r');
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




y_r = zeros(1,n);

% BEC: erasure prob. = 0.05: 0.05 : 0.40;
for max_iterations = 30 : 30 : 120
    num = 0;
    total_trial = 100000;
for erasure_prob = 0.02 : 0.02 : 0.40
    failure = 0;
    num = num + 1;
    for trial = 1 : total_trial
        u = randi([0,1], 1, k);
        x_s = mod(u*G, 2);
        for i = 1 : n
            p = rand;
            if (p <= erasure_prob)
                y_r(i) = 2;
            else 
                y_r(i) = x_s(i);
            end
        end
        
        [decoded_x, convergence] = SPA_BEC(y_r, n, m, row_w, col_w, variable, check, max_iterations, P);
        if ~(convergence && isequal(decoded_x, x_s))
            failure = failure + 1;
        end
    end
    worderr(max_iterations/30, num) = failure/total_trial; 
end
end



marks = {'-rp', '--bo', '--kp', '--g^'};



for max_iterations = 30 : 30 : 120
    semilogy( 0.05 : 0.05 : 0.40, worderr(max_iterations/200,:), marks{max_iterations/200});
end

hold on;

title('Word Error Rate vs erasure probabilty under Binary Erasur Channel using the Sum-Product Algorithm(#Total trial=1000000)');
xlabel('Erasure probabilty');
ylabel('Word Error Rate');

legend('max_iteration = 200', 'max_iteration = 400', 'max_iteration = 600', 'max_iteration = 800');

hold off;