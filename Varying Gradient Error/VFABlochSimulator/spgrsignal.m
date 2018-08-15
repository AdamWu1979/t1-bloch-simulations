%%	function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	following Nex excitations.
%
function [Msig,Mss]=spgrsignal(flip,T1,T2,TE,TR,PartialDephasing,df,Nex,inc)


%% Set undefined function-required variables.
%

if (nargin < 8)
	inc = 117/180*pi;
end;
if (nargin < 7)
	Nex = 100;
end;
if (nargin < 6)
	df = 0;
end;

%% Set up spin properties
%

Nf = 100;	% Simulate 100 different gradient-spoiled spins.
phi = ((1-Nf/2):Nf/2)/Nf*2*pi*PartialDephasing; % Radian phase vector going from 2Pi/Nf to 2Pi in 2Pi/Nf increments. (WHAT DOES IT _REPRESENT_ ?) _MJB

%%**** The definition of this M is obsolete due to the fact that it's
%%overriden  just a few lines over (the third dimension disapears, Nex is
%%no longer a part of this vector) - *** Should be deleted. _MJB
M=zeros(3,Nf,Nex+1); % Creates magnetization vector for each spin, for each excitation.


%% Calculate free-precession matrices
%

%"A" is decay and phase gained due to off resonance, "B" is regrowth
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)]; % Sets initial magnetization for every spin [0;0;1]
on = ones(1,Nf); % Vector to ensure size of matrices in further calculations 
	
Rfph = 0;       % Rf phase
Rfinc = inc;    %*****WRONG****"inc" is in ##deg## here when run by default VFA_simulation.m*******

for n=1:Nex

    % *Not sure why it's necessayr to define A and B here, could have skipped that and went straight to calculating M. _MJB   
	A = Ate * throt(flip,Rfph);
	B = Bte;
	M = A*M+B*on; % M is rotated, then decayed for TE. Regrowth factor is added.

	Msig = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph); % Complex signal by adding up all the spins
	Mss = M; % At ths end of the simulation, this complexe magnetization matrix for all spins is outputed by the function.
    
    % *See, like this! _MJB
	M=Atr*M+Btr*on; % Relaxation during rest of TR after TE

    % What is this? Fake spin dephasing? Why not use frequency bandwidth from T2 to calculate the theoretical spin dephasing that occurs in TR? _MJB
	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
	end;

    % To make sure spoiling is ideal
    % M(1:2, :) = 0; 
    
	Rfph = Rfph+Rfinc; % Calculate the next RF phase
	Rfinc = Rfinc+inc; % Calculate the next RF increment
end;


		







