%% Variabl Flip Angle Bloch Simulator (4 Flips Angles) - Varying The Incomplete Spoiling
% This simulator simulated the signal amplitude from the VFA pulse
% sequence for incomplete spoiling (from gradient spoiling *only*)
%
%
% Main code author: Mathieu Boudreau
% Bloch code author: Mathieu Boudreau, Nikola Stikov
% T1 fitting code author: Nikola Stikov
% Date: November 2012


%% Clear Matlab Session

close all
clear all
clc

%% Parameter initialization
%

% In ms
T1 = 825.5; 
T2 = 100;  
TE = 3.5;  
TR = 15; 

df = 0;    % On resonance
Nex = 200; % Number of excitation to achieve steady state
inc = 117; % Phase spoiling in deg.
color = {'r', 'b', 'm', 'k'};

flipAngle = [3/180*pi 10/180*pi 20/180*pi 30/180*pi];  % FA in radians

B1mapError = 0.95; % Ratio of how the B1 corrected alpha (alpha) differs from the actual implemented pulse.

PartialDephasing = 0.80:0.01:1; % The variable *that is changed* for this simulations.

%% Calculate/Fit signal
%

for kk = 1:length(PartialDephasing)    
    %% Calculate SPGR signal with RF spoiling for each FA
    %
    [Msig1,Mss]=spgrsignal(flipAngle(1),T1,T2,TE,TR,PartialDephasing(kk),df,Nex,inc);
    [Msig2,Mss]=spgrsignal(flipAngle(2),T1,T2,TE,TR,PartialDephasing(kk),df,Nex,inc);
    [Msig3,Mss]=spgrsignal(flipAngle(3),T1,T2,TE,TR,PartialDephasing(kk),df,Nex,inc);
    [Msig4,Mss]=spgrsignal(flipAngle(4),T1,T2,TE,TR,PartialDephasing(kk),df,Nex,inc);
    
    % Signal
    VFA(1) = abs(Msig1);
    VFA(2) = abs(Msig2);
    VFA(3) = abs(Msig3);
    VFA(4) = abs(Msig4);

    %% T1 Fitting
    %
    
    xVFAFitTerm = VFA./tan(flipAngle*B1mapError);
    yVFAFitTerm = VFA./sin(flipAngle*B1mapError);

    fitResults=polyfit(squeeze(xVFAFitTerm), squeeze(yVFAFitTerm), 1); % Linear fit
    slope = fitResults(1); %fitResults(2) is the y-intercept

    % This is a bad name for T1, because RF spoiling is used - should say not perfect spoiling. _MJB    
    t1Unspoiled(kk) = abs(-TR./log(slope)); % Without noise
end


%% Plot Figures
%
figure; 
hold on; 
plot(PartialDephasing, T1*ones(numel(PartialDephasing), 1), 'k', PartialDephasing, t1Unspoiled, 'b');
