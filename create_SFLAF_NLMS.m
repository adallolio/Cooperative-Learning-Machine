function F = create_SFLAF_NLMS(M, muL, muFL, exord, memord, Mi, Me, deltaL, deltaFL, type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE STRUCTURE: SPLIT FUNCTIONAL LINK ADAPTIVE FILTER (SFLAF)
% USING NORMALIZED LEAST MEAN SQUARES (NLMS) ALGORITHM.

% DESCRIPTION:
% Linear combination of a linear NLMS filter and a nonlinear filter.
% Both the filters are updated using the overall error signal.
%
%
% INPUT PARAMETERS:
%   M: filter length
%   muL: step size of the linear filter
%   muFL: step size of the NLMS-FLAF nonlinear filter
%   exord: expansion order
%   memord: memory order - memord = 0 for memoryless functional links
%   Mi: length of the FL input buffer selected for the expansion
%   Me: length of the expanded buffer
%   deltaL: regularization parameter of the linear filter
%   deltaFL: regularization parameter of the NLMS-FLAF nonlinear filter
%   type: functional expansion type
%
% OUTPUT PARAMETERS:
%   F.M: filter length
%   F.muL: step size of the NLMS filter
%   F.muFL: step size of the NLMS-FLAF nonlinear filter
%   F.exord: expansion order
%   F.memord: memory order - memord = 0 for memoryless functional links
%   F.Mi: length of the FL input buffer selected for the expansion
%   F.Me: length of the expanded buffer
%   F.deltaL: regularization parameter of the NLMS filter
%   F.deltaFL: regularization parameter of the FL nonlinear filter
%   F.exord2: number of different sin and cos for each sample
%   F.xBuff: input buffer [M x 1]
%   F.wL: coefficient vector for the NLMS filter [M x 1]
%   F.wFL: coefficient vector for the FL filter [Me x 1]
%   F.exp_type: index of functional expansion type
%
%
% REFERENCES:
%   [1] D. Comminiello, L.A. Azpicueta-Ruiz, M. Scarpiniti, A. Uncini and 
%       J. Arenas- García, "Functional Link Adaptive Filters for 
%       Nonlinear Acoustic Echo Cancellation", in IEEE Transactions on 
%       Audio, Speech and Language Processing, vol. 21, no. 7, pp. 
%       1502-1512, July 2013.
%   [2] D. Comminiello, L.A. Azpicueta-Ruiz, M. Scarpiniti, A. Uncini and 
%       J. Arenas- García, "Functional Links Based Architectures for 
%       Nonlinear Acoustic Echo Cancellation", HSCMA' 11, Edinburh, UK, 
%       May 2011.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 10              % NARGIN = number argument input
    type = 'Tri';                                                          % Trigonometric expansion as default expansion type
end
if strcmp('Tri',type)
    exp_type = 1;
elseif strcmp('Cheb',type)
    exp_type = 2;
elseif strcmp('Leg',type)
    exp_type = 3;
else
    error('Invalid expansion type!')
end
if nargin < 8
    deltaL = 1e-3;                                                         % Default value for the regularization factor of the NLMS filter
end
if nargin < 9
    deltaFL = deltaL;                                                      % Default value for the regularization factor of the NLMS-FLAF filter
end
if nargin < 6
    Mi = round((1/8)*M);                                                   % Default length of the FL input buffer selected for the expansion
end
if nargin < 5
    memord = 0;                                                                                     % Memoryless functional link for default
end
if nargin < 4
    exord = 5;                                                                                      % Default FLAF expansion order
end
if nargin < 3
    muFL = 0.1;                                                                                     % Default step size value for the NLMS-FLAF filter
end
if nargin < 2
    muL = 0.2;                                                                                      % Default step size value for the NLMS filter
end
exord2 = 2*exord;                                                                                   % Number of sin and cos for trigonometric expansion
if nargin < 7
    if exp_type == 1
        Me = Mi*(exord2) + (Mi-memord)*(memord*exord2) + sum(((2:memord)-1)*exord2);                % Length of the expanded buffer for trigonometric expansion
    else
        Me = Mi*exord;                                                                              % Length of the expanded buffer for Chebyshev and Legendre expansion
    end
end
xBuff = zeros(M,1);                                                        % Input buffer inizialization
wL = zeros(M,1);                                                           % Coefficient vector inizialization for the NLMS filter
wFL = zeros(Me,1);                                                         % Coefficient vector inizialization for the NLMS-FLAF nonlinear filter

F = struct('M', M, 'muL', muL, 'muFL', muFL, 'exord', exord,...
    'memord', memord, 'Mi', Mi, 'Me', Me, 'deltaL', deltaL,...
    'deltaFL', deltaFL, 'exord2', exord2, 'xBuff', xBuff, 'wL',...
    wL, 'wFL', wFL, 'exp_type', exp_type);                                                                             % Structure of the NLMS adaptive filter