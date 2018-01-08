function F = create_AF_NLMS(M, mu, delta, wi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE STRUCTURE: NORMALIZED LEAST MEAN SQUARES (NLMS) ADAPTIVE FILTER.

% INPUT PARAMETERS:
%   w0: impulse response
%   M: filter length
%   mu: step size
%   delta: regularization factor
%   wi: initial value for filter coefficients
%
% OUTPUT PARAMETERS:
%   F.w0: impulse response
%   F.M: filter length
%   F.mu: step size
%   F.delta: regularization factor
%   F.xBuff: input buffer vector at time instant n=0 [M,1]
%   F.w: filter vector at time instant n=0 [M,1]
%
% REFERENCES: 
%   [1] A. Uncini, "Elaborazione Adattativa dei Segnali", Aracne Editrice,
%       ISBN: 978:88-548-3142-I, 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 4
    wi = 0;                                                                % Initial value for the filter coefficients
end
xBuff = zeros(M,1);                                                        % Input buffer inizialization
w = wi*ones(M,1);                                                          % Filter vector inizialization

F = struct('M', M, 'mu', mu, 'delta', delta, 'xBuff', xBuff,...
           'w', w);                                                        % Structure of the NLMS adaptive filter