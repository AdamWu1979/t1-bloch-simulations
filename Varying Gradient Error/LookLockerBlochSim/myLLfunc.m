function func = myLLfunc(x,params)
%MYVFAORIGFUNCTION Summary of this function goes here
%   Detailed explanation goes here
    
    alpha = params.alpha;
    beta = params.beta;
    TI1 = params.TI1;
    TI2 = params.TI2;
    TR = params.TR;
    Nll = params.Nll;
    
   calpha = cos(alpha);
   cbeta = cos(beta);
   
   tr = TR - TI1 - (Nll-1).*TI2;
   
   E1 = exp(-TI1./x(2));
   E2 = exp(-TI2./x(2));
   Er = exp(-tr/x(2));
   
   F = (1-E2)./(1-cbeta.*E2);
   Qnom = F.*calpha.*cbeta.*Er.*E1.*(1-(cbeta.*E2).^(Nll-1)) + calpha.*E1.*(1-Er)-E1+1;
   Qdenom = 1-calpha.*cbeta.*Er.*E1.*(cbeta.*E2).^(Nll-1);
   Q=Qnom/Qdenom;
   
   Meq = 1;
   
   % Array pre-loading
   Mz = zeros(Nll,1);
   Mss = zeros(Nll,1);
   
   %Calculate steady state longitudinal magnetization and signal
   for ii=1:Nll

        Mz(ii) = F+(cbeta.*E2).^(ii-1).*(Q-F);
        func(ii) = abs(x(1).*sin(beta).*Mz(ii));
   end
end

