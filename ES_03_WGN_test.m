% Estimation of Correlation matrix of White Gaussian Noise (WGN) 
%
clear all; close all;
N=5;
N_of_Exp = 100000;
R = zeros(N,N);
StdDev = sqrt(1); % sqrt(Sigma^2);
for n = 1 : N_of_Exp
    x = StdDev*randn(N,1);
    R = R + x*x';
end
R./N_of_Exp;
