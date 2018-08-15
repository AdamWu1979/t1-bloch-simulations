%% Inversion Recovery Bloch Simulator - Varying The Flip Angle Error
% This simulator simulated the signal amplitude from the inversion recovery
% pulse sequence for a range of flip angle error (after B1 correction)
%
%
% Main code author: Mathieu Boudreau
% Bloch code author: Mathieu Boudreau, Nikola Stikov
% T1 fitting code author: Joel Barral
% Date: November 2012


%% Clear Matlab Session
%

clear all

%% Code Flags
%

% crusherFlag = 1 -> complete spoiling 
% crusherFlag = 2 -> partial spoiling
crusherFlag = 2; 

%% Parameter initialization
%

load('defaultLL.mat')   % All pulse sequence parameters are set here.

T1est = 900.5;          % Estimate of T1 value when fitting

% Alpha = inversion pulse; Beta = excitation pulse
defaultAlpha = alpha;   % Default here actually the nominal flip angles in radians (before B1 correction)
defaultBeta = beta;     

FAconst = defaultAlpha/defaultBeta; % Constant ratio between all alpha and betas.


alpha = 160;           % The variable that is changed for this simulation. (Value calculated AFTER B1 correction)
beta = alpha/FAconst;

alpha = deg2rad(alpha);
beta = deg2rad(beta);

B1mapError = 0.90:0.01:1; % Ratio of how the B1 corrected alpha (alpha) differs from the actual implemented pulse.

%% ***Run bloch simulator***
%

for jj = 1:length(TI) % Loop over TI acquisitions
    
    % Calculate Look Locker signal with RF spoiling for each FA
    [Msig,Mz]=IRsignal(alpha, beta,TI(jj),T1,T2,TE,TR,crusherFlag,df,Nex,inc); % Alpha changed inside this function
    
    simMss(jj) = Msig; % Signal
    
end

%% Fit for T1
%

% Set fitting parameters
extra.TR = TR;  % Repetition time (TR)
extra.T1Vec = 1:5000;  % Initial grid points for the T1 search
extra.tVec = TI; % Inversion times (TIs) considered
extra.kInit =  2; % Is this even used? What is it?
extra.T1Init =  200; % Initial T1 guess?

nlsS = getNLSStruct(extra);          

data = simMss; % This selects TI = 50, 400, 110, 2500

% Fitting is done here
[T1Hat, bHat, aHat, residual] = rdNls(data, nlsS);
%[T1Hat, bHat, aHat, residual] = lmSphPr(abs(data), extra);

%% Plot figures
%

figure(), plot(B1mapError,T1Hat, 'b'), hold on, plot(B1mapError,ones(length(B1mapError),1)*T1, 'r')
xlabel('Alpha flip angle (deg)')
ylabel('T1 (ms)')
title('Fitted T1 values for varying alpha values WITH CRUSHING, with the fitting procedure assuming the nominal beta value')
legend('Calculated T1','True T1')

