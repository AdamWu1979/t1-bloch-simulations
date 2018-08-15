%%	function [Msig,Mss]=LLsignal(alpha, beta, flip,TI1,TI2,T1,T2,TE,TR,Nll,df,Nex,inc)
%
%	Function calculates the signal from an RF-spoiled sequence
%	following Nex excitations.
%
function [Msig,MLong]=IRsignal(alpha, beta,TI,T1,T2,TE,TR,crushFlag, PartialDephasing, df,Nex,inc)


%% Set undefined function-required variables.
%

if (nargin < 12)
	inc = 0;
end;
if (nargin < 11)
	Nex = 100;
end;
if (nargin < 10)
	df = 0;
end;

%% Set up spin properties
%

Nf = 100;	% Simulate 100 different gradient-spoiled spins.
phi = (1:Nf)/Nf*2*pi; % Radian phase vector going from 2Pi/Nf to 2Pi in 2Pi/Nf increments.

%% Calculate free-precession matrices
%

%"A" is decay and phase gained due to off resonance, "B" is regrowth
[Ate,Bte] = freeprecess(TE,T1,T2,df);                       % Magnetization decayed (A) and regrowth (B) between the beta pulse and measurement.
[Ati,Bti] = freeprecess(TI,T1,T2,df);                    % Magnetization decayed (A) and regrowth (B) between the alpha pulse and beta pulse.
[Atr,Btr] = freeprecess(TR-TI-TE,T1,T2,df);       % Magnetization decayed (A) and regrowth (B) between the LAST measurement and the next TR.

M = [zeros(2,Nf);ones(1,Nf)]; % Sets initial magnetization for every spin [0;0;1]
on = ones(1,Nf); % Vector to ensure size of matrices in further calculations 
	
Rfph = 0;       % Rf phase
Rfinc = inc;    

for n=1:Nex %nth TR
    
    A1 = Ati * throt(alpha, Rfph);
    B1 = Bti;
    
    % Apply alpha pulse, then decay/regrow until beta pulse
    M = A1*M+B1*on;
    
    if crushFlag == 1
       M(1:2,:) = 0; %Crush signal at the end of TI1, but not during the TI2 
    elseif crushFlag == 2
       phi2 = ((1-Nf/2):Nf/2)/Nf*2*pi*PartialDephasing;
       for k=1:Nf
       	   M(:,k) = zrot(phi2(k))*M(:,k); % Dephase spins.
       end
    end

    A2 = Ate * throt(beta, Rfph);
    B2 = Bte;

    % Apply beta pulse, then decay/regrow until measurement time.
    M=A2*M+B2*on; 

    Msig = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph); % Complexe transverse magnetization.
    MLong = mean(M(3,:)); % Longitudinal magnetization at measurement time.

    A3 = Atr;
    B3 = Btr;

    % Decay/regrow  from the last measurement until the start of the next TR.
    M=A3*M+B3*on;
    
    M(1:2,:) = 0;   
    Rfph = Rfph+Rfinc; % Calculate the next RF phase
    Rfinc = Rfinc+inc; % Calculate the next RF increment

end;


		







