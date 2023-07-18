clc;clear

order = 10;
[lambda,weight] = quadpts1(order);

weight = sort(weight);
N = length(weight);

G = repmat(weight,N,1)

G*G'