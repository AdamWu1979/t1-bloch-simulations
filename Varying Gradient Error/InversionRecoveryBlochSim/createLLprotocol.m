clear all

alpha = 180;        % In deg - Variable parameter
alpha = deg2rad(alpha); % Changed into radians

beta = 90;
beta = deg2rad(beta); % Changed into radians

% In ms
T1 = 825.5;          
T2 = 100;           
TE = 11;

TI = [30 530 1030 1530];
TR = 1550;           

df = 0;             % Off-resonance frequency

Nex = 200;          % # to induce SS
inc = 0; %120; %117; % Spoiling phase inc
Nll = 4;

save('defaultLL.mat')