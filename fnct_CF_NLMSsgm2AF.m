function [e, y, F] = fnct_CF_NLMSsgm2AF(x, d, F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON-LINE FUNCTION: ADAPTIVE FILTER COMBINING TWO ADAPTIVE FILTER OUTPUTS
% BY USING THE SIGMOID FUNCTION.

% INPUT PARAMETERS:
%   x: input vector containing the AF outputs at n-th time instant
%   d: desired signal sample at n-th time instant
%   F.K: combining filter length
%   F.mu: step size
%   F.beta: smoothing factor for power normalization
%   F.h: combining filter coefficient at time instant n-1
%   F.a: auxiliary adaptation coefficient at time instant n-1
%   F.r: power normalization factor at time instant n-1

%
% OUTPUT PARAMETERS:
%   e: error signal sample at n-th time instant
%   y: output signal sample at n-th time instant
%   F.K: combining filter length
%   F.mu: step size
%   F.beta: smoothing factor for power normalization
%   F.h: combining filter coefficient at n-th time instant
%   F.a: auxiliary adaptation coefficient at n-th time instant
%   F.r: power normalization factor at n-th time instant
%
%
% REFERENCES: 
%   [1] L.A. Azpicueta-Ruiz, A.R. Figueiras-Vidal, and J. Arenas-García, "A
%       Normalized Adaptation Scheme for the Convex Combination of Two
%       Adaptive Filters", in Proc. of the IEEE Int. Conf. Acoust. Speech
%       Signal Process. (ICASSP '08), pp. 3301-3304, 2008.
%   [2] A. Uncini, "Elaborazione Adattativa dei Segnali", Aracne Editrice,
%       ISBN: 978:88-548-3142-I, 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------Compute the mixing parameter---------------------------------------------%
if abs(F.a) > 4                                                                     % Limit the adaptation parameter
    F.a = sign(F.a)*4;
end
F.h = logsig(F.a);                                                                  % Compute the sigmoid function
F.r = F.beta*F.r + (1 - F.beta)*(x(1) - x(2))^2;                                    % Power normalization estimation

%----------Compute overall output and error signals---------------------------------%
y = F.h*x(1) + (1 - F.h)*x(2);                                                      % Overall output signal
e = d - y;                                                                          % Overall error signal

%----------Update the adaptation parameter------------------------------------------%
F.a = F.a + (F.mu/F.r)*F.h*(1 - F.h)*e*(x(1) - x(2));                               % NLMS adaptation of the mixing parameter