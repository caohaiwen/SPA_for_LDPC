function ReadOut_Plot_BEC(data)
% read out the BER & WER data for BEC


file = fopen(data, 'r');
num_line = 0;
num = 0;
its = 0;
prv_its = 0;
erasure_prob = 0.02 : 0.02 :0.52;

current = 0;

while ~feof(file)
    num_line = num_line + 1;
    current = current + 1;
    content = fgetl(file);
    s = str2num(content);
    
    its = s(1);
    
    if (prv_its ~= its)
        current = 1;
        num = num + 1;
        prv_its = its;
        iteration(num) = s(1);
    end
    
    if mod(num_line,2) == 1
       word_err(num,(current + 1)/2) = s(3);
    else 
       bit_err(num, current/2) = s(3);
    end
   
end

fclose(file);

marks = {'-rp', '--bo', '--k^', '-gd', ':mh'};
subplot(2,1,1);


set(gca, 'YScale', 'log');
xticks(erasure_prob);

hold on
for i = 1 : 5
    semilogy(erasure_prob, word_err(i,:), marks{i});
end

title('Word Error Rate vs Erasure Prob under BEC using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('Erasure Probabilty');
ylabel('Word Error Rate');


legend('max_iterations = 8', 'max_iterations = 16', 'max_iterations = 32' , 'max_iterations = 64' , 'max_iterations = 128');
hold off;

subplot(2,1,2);

set(gca, 'YScale', 'log');
xticks(erasure_prob);

hold on

for i = 1 : 5
    semilogy(erasure_prob, bit_err(i,:), marks{i});
end

title('Bit Error Rate vs Erasure Prob under BEC using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('Erasure Probabilty');
ylabel('Bit Error Rate');


legend('max_iterations = 8', 'max_iterations = 16', 'max_iterations = 32' , 'max_iterations = 64' , 'max_iterations = 128');

hold off;


