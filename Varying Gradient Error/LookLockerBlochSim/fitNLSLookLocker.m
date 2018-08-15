function [fittedT1,fittedConst]= fitNLSLookLocker(simMss,alpha, beta,TI1,TI2,T1est,TR,Nll)
%FITNLSVFA Summary of this function goes here
%   Detailed explanation goes here

    params.alpha = alpha;
    params.beta = beta;
    params.TI1 = TI1;
    params.TI2 = TI2;
    params.TR = TR;
    params.Nll = Nll;
    
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    xdata = [30, 530, 1030, 1530];
    fitParamLL=lsqcurvefit(@(x,params)myLLfunc(x,params), [0.5, T1est], params, simMss',[],[],opts);

    
    
    % m = exp(-TR/T1), solve for T1
    fittedT1 = fitParamLL(:,2);

    % b = M0(1-exp(-TR/T1)), solve for M0
    fittedConst = fitParamLL(:,1);

end

