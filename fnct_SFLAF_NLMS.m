function [e, y, yFL, F] = fnct_SFLAF_NLMS(x, d, F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON-LINE FUNCTION: SPLIT FUNCTIONAL LINK ADAPTIVE FILTER (SFLAF)
% USING NORMALIZED LEAST MEAN SQUARES (NLMS) ALGORITHM.

% DESCRIPTION:
% Linear combination of a linear NLMS filter and a nonlinear filter.
% Both the filters are updated using the overall error signal.
%
%
% INPUT PARAMETERS:
%   x: input signal sample at n-th time instant
%   d: desired signal sample at n-th time instant
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
% OUTPUT PARAMETERS:
%   e: error signal sample at n-th time instant
%   y: output signal sample at n-th time instant
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


%----------Update input signal buffer-----------------------------------------------%
F.xBuff(2:F.M) = F.xBuff(1:F.M-1);                                                  % Shift temporary input signal buffer down
F.xBuff(1) = x;                                                                     % Assign current input signal sample

%----------Functional link expansion------------------------------------------------%
g = flex(F.xBuff, F.Mi, F.Me, F.exord, F.memord, F.exp_type);                       % Generate the expanded buffer

%----------Compute individual output signals----------------------------------------%
normL = F.xBuff'*F.xBuff + F.deltaL;                                                % Norm of the input layer vector for the NLMS filter
normFL = g'*g + F.deltaFL;                                                          % Norm of the input layer vector for the NLMS-FLAF filter
yL = F.wL'*F.xBuff;                                                                 % Compute the output signal for the NLMS filter
yFL = F.wFL'*g;                                                                     % Compute the output signal for the NLMS-FLAF filter

%----------Compute overall output and error signals---------------------------------%
y = yL + yFL;                                                                       % Combination of individual output signals
e = d - y;                                                                          % Overall error signal

%----------Update the filters-------------------------------------------------------%
F.wL = F.wL + (F.muL/normL)*e*F.xBuff;                                              % Update of the NLMS filter
F.wFL = F.wFL + (F.muFL/normFL)*e*g;                                                % Update of the NLMS-FLAF filter

