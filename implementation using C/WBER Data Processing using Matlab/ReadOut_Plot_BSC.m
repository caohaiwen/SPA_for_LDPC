% read out the BER & WER data for AWGN

clc;
clear;
file = fopen("./PCMatrix(96.3.963 (N=96,K=48,M=48,R=0.5))/WBER under BSC.txt", 'r');
num_line = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       cross_prob((num_line + 1)/2) = s(1);
       word_err(1,(num_line + 1)/2) = s(2);
   else 
       bit_err(1, num_line/2) = s(2);
   end
   
end

fclose(file);

len = length(cross_prob);

file = fopen("./PCMatrix(N=204,K=102,M=102,R=0.5)/WBER under BSC.txt", 'r');
num_line = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       word_err(2,(num_line + 1)/2) = s(2);
   else 
       bit_err(2, num_line/2) = s(2);
   end
   
end

word_err(2, :) = word_err(2, :) + [zeros(1, num_line /2), ones(1, len - num_line / 2)];
bit_err(2, :) = bit_err(2, :) + [zeros(1, num_line /2), ones(1, len - num_line / 2)] / 2;

fclose(file);

file = fopen("./PCMatrix(816.3.174 (N=816,K=408,M=408,R=0.5))/WBER under BSC.txt", 'r');
num_line = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       word_err(3,(num_line + 1)/2) = s(2);
   else 
       bit_err(3, num_line/2) = s(2);
   end
   
end


word_err(3, :) = word_err(3, :) + [zeros(1, num_line /2), ones(1, len - num_line / 2)];
bit_err(3, :) = bit_err(3, :) + [zeros(1, num_line /2), ones(1, len - num_line / 2)] / 2;

fclose(file);

marks = {'-rp', '--bo', '--kp'};
subplot(2,1,1);

set(gca, 'YScale', 'log');
xticks(cross_prob);

hold on
for i = 1 : 3
    semilogy(cross_prob, word_err(i,:), marks{i});
end

title('Word Error Rate vs Crossover Prob. under BSC using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('Crossover Probablity');
ylabel('Word Error Rate');


legend('N=96,K=48,M=48,R=0.5', 'N=204,K=102,M=102,R=0.5', 'N=816,K=408,M=408,R=0.5');

hold off;

subplot(2,1,2);

set(gca, 'YScale', 'log');
xticks(cross_prob);

hold on
for i = 1 : 3
    semilogy(cross_prob, bit_err(i, :), marks{i});
end


title('Bit Error Rate vs Crossover Prob. under BSC using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('Crossover Probablity');
ylabel('Bit Error Rate');

legend('N=96,K=48,M=48,R=0.5', 'N=204,K=102,M=102,R=0.5', 'N=816,K=408,M=408,R=0.5');

hold off;

