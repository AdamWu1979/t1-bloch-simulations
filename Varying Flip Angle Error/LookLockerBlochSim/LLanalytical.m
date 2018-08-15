function [Mss, Mz]=LLanalytical(alpha, beta,TI1,TI2,T1,T2,TE,TR,Nll)
%LLANALYTICAL Calculates the signal and 
    
   calpha = cos(alpha);
   cbeta = cos(beta);
   
   tr = TR - TI1 - (Nll-1).*TI2;
   
   E1 = exp(-TI1./T1); 
   E2 = exp(-TI2./T1);
   Er = exp(-tr/T1);
   
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
        Mz(ii) = Meq.*(F+(cbeta.*E2).^(ii-1).*(Q-F));
        Mss(ii) = abs(sin(beta).*Mz(ii).*exp(-TE/T2));
   end
end

