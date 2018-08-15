%% Look-Locker Bloch Simulator - Varying The Incomplete Spoiling
% This simulator simulated the signal amplitude from the Look-Locker
% pulse sequence for incomplete spoiling (from gradient spoiling *only*)
%
%
% Main code author: Mathieu Boudreau
% Bloch code author: Mathieu Boudreau, Nikola Stikov
% T1 fitting code author: Mathieu Boudreau
% Date: November 2012


%% Clear Matlab Session
%

clear all;
close all;
clc

%% Code Flags
%

% crusherFlag = 1 -> complete spoiling 
% crusherFlag = 2 -> partial spoiling
crusherFlag = 2; 

%% Parameter initialization
%

load('LLprotocol.mat')

T1est = 900.5;          % Estimate of T1 value when fitting

% Alpha = inversion pulse; Beta = excitation pulse
defaultAlpha = alpha;   % Default here actually the nominal flip angles in radians (before B1 correction)
defaultBeta = beta;     

FAconst = defaultAlpha/defaultBeta; % Constant ratio between all alpha and betas.


alpha = 160;           % The variable that is changed for this simulation. (Value calculated AFTER B1 correction)
beta = alpha/FAconst;

alpha = deg2rad(alpha);
beta = deg2rad(beta);

B1mapError = 0.95; % Ratio of how the B1 corrected alpha (alpha) differs from the actual implemented pulse.

PartialDephasing = 0.80:0.01:1; % The variable *that is changed* for this simulations.

%% ***Run bloch simulator***
%

for jj=1:length(PartialDephasing) % Loop over spoiling error
    
    % Calculate Look Locker signal with RF spoiling for each FA
    [Msig,Mz]=LLsignal(alpha, beta,TI1,TI2,T1,T2,TE,TR,crusherFlag,PartialDephasing(jj),Nll,df,Nex,inc);
    simMss(:,jj) = abs(Msig); % Signal
    simMz(:,jj) = Mz;         % Longitudinal magnetization
    
end

%% Fit for T1
%

for kk = 1:length(PartialDephasing) % Loop over Spoiling error
    
    [fittedT1,fittedConst]= fitNLSLookLocker(simMss(:,kk),alpha*B1mapError, beta*B1mapError,TI1,TI2,T1est,TR,Nll);
    t1Values(kk) = fittedT1;
    
end

%% Plot Figures
%

figure(), plot(PartialDephasing,t1Values, 'b'), hold on, plot(PartialDephasing,ones(length(PartialDephasing),1)*T1, 'r')
xlabel('Alpha flip angle (deg)')
ylabel('T1 (ms)')
title('Fitted T1 values for varying alpha values WITH CRUSHING, with the fitting procedure assuming the nominal beta value')
legend('Calculated T1','True T1')