function F = create_CF_NLMSsgm2AF(K, mu, wi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE STRUCTURE: ADAPTIVE FILTER COMBINING TWO ADAPTIVE FILTER OUTPUTS
% BY USING THE SIGMOID FUNCTION.

% INPUT PARAMETERS:
%   K: combining filter length
%   mu: step size
%   delta: regularization factor
%   wi: initial value for filter coefficients
%
% OUTPUT PARAMETERS:
%   F.K: combining filter length
%   F.mu: step size
%   F.beta: smoothing factor for power normalization
%   F.h: combining filter coefficients at time instant n=0
%   F.a: auxiliary adaptation coefficients at time instant n=0
%   F.r: power normalization factor at time instant n=0
%
% REFERENCES: 
%   [1] L.A. Azpicueta-Ruiz, A.R. Figueiras-Vidal, and J. Arenas-García, "A
%       Normalized Adaptation Scheme for the Convex Combination of Two
%       Adaptive Filters", in Proc. of the IEEE Int. Conf. Acoust. Speech
%       Signal Process. (ICASSP '08), pp. 3301-3304, 2008.
%   [2] A. Uncini, "Elaborazione Adattativa dei Segnali", Aracne Editrice,
%       ISBN: 978:88-548-3142-I, 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    wi = 0.5;                                                              % Initial value for the filter coefficients
end
if nargin < 2
    mu = 0.5;                                                              % Default step size parameter
end

h = wi;                                                                    % Combining filter vector inizialization
a = 4;                                                                     % Adaptation parameter initialization
r = 1;                                                                     % Power factor initialization
beta = 0.9;                                                                % Power combination factor

F = struct('K', K, 'mu', mu, 'beta', beta, 'r', r, 'a', a, 'h', h);        % Structure of the NLMS adaptive filter