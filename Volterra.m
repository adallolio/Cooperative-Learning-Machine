%%   Volterra_LMS.m
%       Implements the Volterra LMS algorithm for REAL valued data.
%       (Algorithm 11.1 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%  
%   Input parameters:
%       N      : members of ensemble.
%       it     : iterations.
%       Sx     : standard deviation of input.
%       Sn     : standard deviation of measurement noise.
%       u      : convergence factor matrix.
%       Nw     : length of the adaptive filter.
%      
%   Output parameters:
%       MSE    : mean-square error.
%%


clear all;		


% Input: 
N      = 1000;
it     = 5e2;
Sx     = 1; 
Sn     = 1e-1;
Nw     = 9;
u      = diag([0.08*ones(1,3) 0.01*ones(1,Nw-3)]);


for j=1:N
   n=Sn*randn(it,1);        % noise at system output 
   x=Sx*randn(it,1);        % input signal
   xl1=zeros(it,1); xl2=xl1; 
   xl1(2:it)=x(1:it-1);    % x(k-1)
   xl2(3:it)=x(1:it-2);    % x(k-2)
   d=zeros(it,1);
   d=-.76*x-xl1+xl2+.5*x.^2+2*x.*xl2-1.6*xl1.^2+1.2*xl2.^2+.8*xl1.*xl2+n; ...
     % unknown system output
   w=zeros(Nw,it);           % initial coefficient vector
   uxl=[x xl1 xl2 x.^2 x.*xl1 x.*xl2 xl1.^2 xl1.*xl2 xl2.^2]'; % input vectors
   
   for i=1:it
      e(i)=d(i)-w(:,i)'*uxl(:,i)+n(i);    % error sample
      y(i)=w(:,i)'*uxl(:,i);              % output sample
      w(:,i+1)=w(:,i)+2*u*e(i)*uxl(:,i);  % new coefficient vector
   end
   mse(j,:)=e.^2;
end 

MSE=mean(mse);


% Output:
figure,
plot(10*log10(MSE));
grid;
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');