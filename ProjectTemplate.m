%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COOPERATIVE FUNCTIONAL LINK ADAPTIVE FILTERS FOR NONLINEAR IDENTIFICATION
% Template.

clear all; close all; clc
disp('***************************************************************************************');
disp('***    COOPERATIVE FUNCTIONAL LINK ADAPTIVE FILTERS FOR NONLINEAR IDENTIFICATION    ***');
disp('***                                EXPERIMENT #0                                    ***');
disp('***************************************************************************************');

% In this example a cooperative architecture involving 2 SFLAFs with
% different expansion orders is implemented.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           INITIALIZATION STAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------Initialization of global parameters----------------------------------------------%
Lx = 20000;                                                                                 % Signal length - Use larger length if you need
nRun = 3;                                                                                   % Number of runs for the avarage - Use few numbers for initial simulations and at least 1000 runs for final results (without any smoothing filtering)
Mc = 2;                                                                                     % Number of constituents to combine = Combining filter length
% Write other global parameters if any

%----------Generate the unknown system to identify------------------------------------------%
model = 5;                                                                                  % Choose the unknown system - see sigen.m
M = 7;                                                                                      % Length of the impulse response to identify = Filter length
nf = 1;                                                                                     % Number of fragments (i.e. changes of the impuse response)
% Write other parameters related to the unknown system if any

%----------Parameter initialization for the input signal------------------------------------%
insig = 2;                                                                                  % Choose the input signal: "1" white Gaussian noise, "2" coloured Gaussian noise
W0 = 2*rand(M, nf) - 1;                                                                     % Generate nf impulse responses whose coefficients are comprised in [-1,1]
if insig ==2    % if coloured Gaussian noise
    a = 0.8;                                                                                % Coefficient of the autoregressive model
    B = sqrt(1 - a^2);
    A = [1 -a];
end

%----------Initialization of the linear filter----------------------------------------------%
wiAF = 0;                                                                                   % Initialization of filter coefficients
muL = 0.2;                                                                                  % Default step size parameter for the linear filter
deltaL = 1e-2;                                                                              % Default regularization factor for the linear APA filter
eL = zeros(Lx,1);                                                                           % Error signal initialization for the linear APA filter
yL = zeros(Lx,1);                                                                           % Output signal initialization for the linear APA filter
navemseL = zeros(Lx,1);                                                                     % Non-averaged EMSE buffer initialization for the linear filter

%----------Initialization of the SFLAF filter #1--------------------------------------------%
muNL = 0.1;                                                                                 % Default step size parameter for all FLAF-based filters
deltaNL = 1e-1;                                                                             % Default regularization factor for NLMS filters for all FLAF-based filters
type = 'Tri';                                                                               % Trigonometric expansion
exordSA1 = 5;                                                                               % Default expansion order for all FLAF-based filters
memord = 0;                                                                                 % Default memory order for all FLAF-based filters - Try also with memord > 0
Mi = M;                                                                                     % Default input length for the expanded buffer for all FLAF-based filters - Try also with Mi < M
MeSA1 = Mi*(2*exordSA1) + (Mi-memord)*(memord*2*exordSA1) + sum(((2:memord)-1)*2*exordSA1); % Default output length for the expanded buffer for SFLAF #1
eSA1 = zeros(Lx,1);                                                                         % SFLAF error signal initialization with exordSA1
ySA1 = zeros(Lx,1);                                                                         % SFLAF output signal initialization with exordSA1
navemseSA1 = zeros(Lx,1);                                                                   % Non-averaged EMSE buffer initialization for the SFLAF #1 with exordSA1

%----------Initialization of the SFLAF filter #2--------------------------------------------%
exordSA2 = 10;                                                                              % Default expansion order for all FLAF-based filters
MeSA2 = Mi*(2*exordSA2) + (Mi-memord)*(memord*2*exordSA2) + sum(((2:memord)-1)*2*exordSA2); % Default output length for the expanded buffer for SFLAF #2
eSA2 = zeros(Lx,1);                                                                         % SFLAF error signal initialization with exordSA2
ySA2 = zeros(Lx,1);                                                                         % SFLAF output signal initialization with exordSA2
navemseSA2 = zeros(Lx,1);                                                                   % Non-averaged EMSE buffer initialization for the SFLAF #2 with exordSA2

%----------Initialization of the CooSFLAF filter using sigmoid-based convex NLMS------------%
wiCF = 0.5;                                                                                 % Initialization of filter coefficients
mua = 0.5;                                                                                  % Default step size value for the adaptation parameter
eCCNs = zeros(Lx,1);                                                                        % Error signal initialization for the CooSFLAF filter using convex constrained NLMS with sigmoid function
yCCNs = zeros(Lx,1);                                                                        % Output signal initialization for the CooSFLAF filter using convex constrained NLMS with sigmoid function
hCCNs = zeros(Lx,1);                                                                        % Mixing parameters of the combining filter using convex constrained NLMS with sigmoid function
navemseCCNs = zeros(Lx,1);                                                                  % EMSE auxiliary buffer initialization for the combining filter using convex constrained NLMS with sigmoid function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           MAIN STAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for nR = 1:nRun

    if  (mod(nR,100)==0) 
       fprintf('\n%d / %d', nR, nRun); 
    end

    %----------Generation of input and desired signals--------------------------------------%
    dx = zeros(Lx, 1);                                                                      % Desired signal initialization
    u = wgn(Lx, 1, 0);                                                                      % Generate a white Gaussian noise
    if insig == 1
        x = u;                                                                              % Assign white Gaussian noise
    elseif insig == 2
        x = filter(B, A, u);                                                                % Generate AR1 signal
    end
    [d, e0] = sigen(x, W0, model);                                                          % Generate desired signal and additive noise

    %----------Create structures------------------------------------------------------------%
    F_linear = create_AF_NLMS(M, muL, deltaL, wiAF);                                        % Create the linear filter using the NLMS algorithm
    F_SA1 = create_SFLAF_NLMS(M, muL, muNL, exordSA1, memord, Mi, MeSA1, deltaL, deltaNL, type);    % Create SFLAF adaptive filter #1 with exordSA1
    F_SA2 = create_SFLAF_NLMS(M, muL, muNL, exordSA2, memord, Mi, MeSA2, deltaL, deltaNL, type);    % Create SFLAF adaptive filter #2 with exordSA2
    F_CCNs = create_CF_NLMSsgm2AF(Mc, mua, wiCF);                                           % Create the combining filter using sigmoid-based convex NLMS
    
    %----------Main process-----------------------------------------------------------------%
    for n = 1:Lx

        %----------Perform the linear filtering---------------------------------------------%
        [eL(n), yL(n), F_linear] = fnct_AF_NLMS(x(n), d(n), F_linear);                      % Provide current error and output samples for the linear filter
        navemseL(n) = navemseL(n) + (eL(n) - e0(n))^2;                                      % Compute current EMSE for the linear filter
        
        %----------Perform the SFLAF filtering with exordSA1--------------------------------%
        [eSA1(n), ySA1(n), yFL1, F_SA1] = fnct_SFLAF_NLMS(x(n), d(n), F_SA1);               % Provide current error and output samples for the SFLAF #1 with exordSA1
        navemseSA1(n) = navemseSA1(n) + (eSA1(n) - e0(n))^2;                                % Compute current EMSE for the SFLAF #1 with exordSA1

        %----------Perform the SFLAF filtering with exordSA2--------------------------------%
        [eSA2(n), ySA2(n), yFL2, F_SA2] = fnct_SFLAF_NLMS(x(n), d(n), F_SA2);               % Provide current error and output samples for the SFLAF #2 with exordSA2
        navemseSA2(n) = navemseSA2(n) + (eSA2(n) - e0(n))^2;                                % Compute current EMSE for the SFLAF #2 with exordSA2
        
        %----------Perform the Combining filtering------------------------------------------%
        xc = [ySA1(n) ySA2(n)];                                                             % Input vector to the output stage
        [eCCNs(n), yCCNs(n), F_CCNs] = fnct_CF_NLMSsgm2AF(xc, d(n), F_CCNs);                % Provide current error and output samples for the Cooperative architecture using sigmoid-based convex NLMS
        hCCNs(n) = hCCNs(n) + F_CCNs.h;                                                     % Mixing parameter evolution for the combining filter using sigmoid-based convex NLMS
        navemseCCNs(n) = navemseCCNs(n) + (eCCNs(n) - e0(n))^2;                             % Compute current EMSE for the CooSFLAF-APA filter using sigmoid-based convex NLMS

    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PERFORMANCE MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------Excess Mean Square Error---------------------------------------------------------%
emseL = 10*log10(navemseL/nRun);                                                            % EMSE of the linear APA filter
emseSA1 = 10*log10(navemseSA1/nRun);                                                        % EMSE of the SFLAF-APA filter with KNL1
emseSA2 = 10*log10(navemseSA2/nRun);                                                        % EMSE of the SFLAF-APA filter with KNL2
emseCCNs = 10*log10(navemseCCNs/nRun);                                                      % EMSE of the CooSFLAF-APA filter using sigmoid-based convex NLMS

[bButt,aButt] = butter(4, 0.005);             % 0.005 = cut frequency; 4 = filter order     % Butterworth 2nd order smooth filter - Comment for very large nRun
semseL = filtfilt(bButt, aButt, emseL);                                                     % Smoothed EMSE of the linear filter - Comment for very large nRun
semseSA1 = filtfilt(bButt, aButt, emseSA1);                                                 % Smoothed EMSE of the SFLAF #1 with exordSA1 - Comment for very large nRun
semseSA2 = filtfilt(bButt, aButt, emseSA2);                                                 % Smoothed EMSE of the SFLAF #2 with exordSA2 - Comment for very large nRun
semseCCNs = filtfilt(bButt, aButt, emseCCNs);                                               % Smoothed EMSE of the combining filter using sigmoid-based convex NLMS - Comment for very large nRun

%----------Mixing Parameters Average--------------------------------------------------------%
ahCCNs = hCCNs/nRun;                                                                        % Mixing parameter average of the combining filter using sigmoid-based convex NLMS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try all other plots you consider necessary according to the followings:

figure1 = figure('PaperSize', [20.98 29.68], 'Color', [1 1 1]);
axes1 = axes('Parent', figure1, 'YGrid', 'on', 'XGrid', 'on', 'FontSize', 20);
box(axes1, 'on');
hold(axes1, 'all');
% plot(emseL, 'Parent', axes1, 'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Linear Filter', 'YDataSource', 'emseL');  % Use for very large nRun
plot(semseL, 'Parent', axes1, 'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Linear Filter', 'YDataSource', 'emseL'); % Use for small nRun
% plot(emseSA1, 'LineWidth', 2, 'DisplayName', 'SFLAF #1 with exordSA1', 'YDataSource', 'emseSA1');    % Use for very large nRun
plot(semseSA1, 'LineWidth', 2, 'DisplayName', 'SFLAF #1 with exordSA1', 'YDataSource', 'semseSA1');    % Use for small nRun
% plot(emseSA2, 'LineWidth', 2, 'DisplayName', 'SFLAF #2 with exordSA2', 'YDataSource', 'emseSA2');    % Use for very large nRun
plot(semseSA2, 'LineWidth', 2, 'DisplayName', 'SFLAF #2 with exordSA2', 'YDataSource', 'semseSA2');    % Use for small nRun
% plot(emseCCNs, 'LineWidth', 2, 'DisplayName', 'Cooperative LM', 'YDataSource', 'emseCCNs');   % Use for very large nRun
plot(semseCCNs, 'LineWidth', 2, 'DisplayName', 'Cooperative LM', 'YDataSource', 'emseCCNs');    % Use for small nRun
xlabel('Samples');
ylabel('EMSE [dB]');
legend(axes1,'show');

figure2 = figure('PaperSize', [20.98 29.68], 'Color', [1 1 1]);
axes2 = axes('Parent', figure2, 'YGrid', 'on', 'XGrid', 'on', 'FontSize', 20);
box(axes2, 'on');
hold(axes2, 'all');
plot(ahCCNs, 'LineWidth', 2, 'DisplayName', 'Cooperative LM', 'YDataSource', 'hCCNs');
xlabel('Samples');
ylabel('Mixing parameter');
legend(axes2,'show');
