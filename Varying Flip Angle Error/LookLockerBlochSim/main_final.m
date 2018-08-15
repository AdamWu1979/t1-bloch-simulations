%% Look-Locker Bloch Simulator - Varying The Flip Angle Error
% This simulator simulated the signal amplitude from the Look-Locker
% pulse sequence for a range of flip angle error (after B1 correction)
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


%% ***Run bloch simulator***
%

% Calculate Look Locker signal with RF spoiling for each FA (Looping done
% inside LLsignal
[Msig,Mz] = LLsignal(alpha, beta,TI1,TI2,T1,T2,TE,TR,crusherFlag,Nll,df,Nex,inc);
simMss = abs(Msig); % Signal
simMz = Mz;         % Longitudinal magnetization


%% Fit for T1
%

for jj = 1:length(B1mapError) % Loop over all B1 errors
    [fittedT1,fittedConst]= fitNLSLookLocker(simMss',alpha*B1mapError(jj), beta*B1mapError(jj),TI1,TI2,T1est,TR,Nll);
    t1Values(jj) = fittedT1;
end


%% Plot Figures
%

figure(), plot(B1mapError,t1Values, 'b'), hold on, plot(B1mapError,ones(length(B1mapError),1)*T1, 'r')
xlabel('Alpha flip angle (deg)')
ylabel('T1 (ms)')
title('Fitted T1 values for varying alpha values WITH CRUSHING, with the fitting procedure assuming the nominal beta value')
legend('Calculated T1','True T1')