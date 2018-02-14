function [ entropy ] = Compute_Entropy( prob )
%this function is to compute the Entropy (in bit) of the corresponding probabilty.
% input : the probabilty "prob" you want to compute
% output : the correspoding Entropy (in bit)

if sum(prob) ~= 1
    disp("Input Error: It's not a probability distribution");
    return;
end

n = length(prob);

entropy = 0;
for ii = 1 : n
    entropy = entropy + prob(ii) * log2(prob(ii));
end

entropy = - entropy ;

end

