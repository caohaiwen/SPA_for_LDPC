% read out the BER & WER data for AWGN

clc;
clear;
file = fopen("./PCMatrix(96.3.963 (N=96,K=48,M=48,R=0.5))/WBER under AWGN Channel.txt", 'r');
num_line = 0;
num_row = 0;
num_col = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       SNR_b((num_line + 1)/2) = s(1);
       word_err(1,(num_line + 1)/2) = s(2);
   else 
       bit_err(1, num_line/2) = s(2);
   end
   
end

fclose(file);

file = fopen("./PCMatrix(N=204,K=102,M=102,R=0.5)/WBER under AWGN Channel.txt", 'r');
num_line = 0;
num_row = 0;
num_col = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       SNR_b((num_line + 1)/2) = s(1);
       word_err(2,(num_line + 1)/2) = s(2);
   else 
       bit_err(2, num_line/2) = s(2);
   end
   
end

fclose(file);

file = fopen("./PCMatrix(816.3.174 (N=816,K=408,M=408,R=0.5))/WBER under AWGN Channel.txt", 'r');
num_line = 0;
num_row = 0;
num_col = 0;
while ~feof(file)
    num_line = num_line + 1;
    content = fgetl(file);
    s = str2num(content);
   if mod(num_line,2) == 1
       SNR_b((num_line + 1)/2) = s(1);
       word_err(3,(num_line + 1)/2) = s(2);
   else 
       bit_err(3, num_line/2) = s(2);
   end
   
end

fclose(file);

marks = {'-rp', '--bo', '--kp'};
subplot(2,1,1);

set(gca, 'YScale', 'log');
xticks(SNR_b);

hold on
for i = 1 : 3
    semilogy(SNR_b, word_err(i,:), marks{i});
end

title('Word Error Rate vs SNR_b under AWGN Channel using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('SNR per information bit (dB)');
ylabel('Word Error Rate');


legend('N=96,K=48,M=48,R=0.5', 'N=204,K=102,M=102,R=0.5', 'N=816,K=408,M=408,R=0.5');

hold off;

subplot(2,1,2);

set(gca, 'YScale', 'log');
xticks(SNR_b);

hold on
for i = 1 : 3
    semilogy(SNR_b, bit_err(i, :), marks{i});
end

title('Bit Error Rate vs SNR_b under AWGN Channel using the Sum-Product Algorithm(#Total trial=100000)');
xlabel('SNR per information bit (dB)');
ylabel('Bit Error Rate');

legend('N=96,K=48,M=48,R=0.5', 'N=204,K=102,M=102,R=0.5', 'N=816,K=408,M=408,R=0.5');

hold off;

