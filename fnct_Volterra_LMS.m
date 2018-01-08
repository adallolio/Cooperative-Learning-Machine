function [e, y, F] = fcnt_Volterra_LMS (x, d, F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON-LINE FUNCTION: SPLIT FUNCTIONAL LINK ADAPTIVE FILTER (SFLAF)
% USING NORMALIZED LEAST MEAN SQUARES (NLMS) ALGORITHM.



% INPUT PARAMETERS:
%   x: input signal sample at n-th time instant
%   d: desired signal sample at n-th time instant
%   F.Mi: input buffer length
%   F.Me: expanded buffer length
%   F.mu: step size of the filter
%   F.exord: expansion order
%   F.delta: regularization parameter of the LMS filter
%   F.xBuff: input buffer [M x 1]
%   F.w: coefficient vector for the filter [M x 1]

% OUTPUT PARAMETERS:
%   e: error signal sample at n-th time instant
%   y: output signal sample at n-th time instant
%   F.N: input buffer length
%   F.mu: step size of the filter
%   F.exord: expansion order
%   F.delta: regularization parameter of the LMS filter
%   F.xBuff: input buffer [M x 1]
%   F.w: coefficient vector for the filter [M x 1]
%   F.K: coefficients number

%mudiag = zeros(length(F.w));
%var=0;
%for i = 0:(length(F.K)-1)
%    mudiag[var:var+F.K(i)] = F.mu(i);
%    var=var+F.K(i);
%end

%u = diag(mudiag);


%----------Update input signal buffer-----------------------------------------------%
    F.xBuff(2:F.Mi) = F.xBuff(1:F.Mi-1);                                                % Shift temporary input signal buffer down
    F.xBuff(1) = x;                                                                     % Assign current input signal sample


%----------Volterra expansion------------------------------------------------%
g = flex(F.xBuff, F.Mi, F.Me, F.exord, F.memord, F.exp_type);                       % Generate the expanded buffer




%----------Compute individual output signals----------------------------------------%
e = d - '*F.w

end

