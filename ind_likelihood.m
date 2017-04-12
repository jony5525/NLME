%Calculates the individual ln likelihood
function li = ind_likelihood(model,errormodel,tdata,cdata,theta,sigma,eta,bInter,lC,det_res_var,bUDDLike)

if (~bUDDLike)
    ipred = model(tdata,cdata,theta,eta);
    res = tdata(:,3)-ipred; %Individual residuals

    if (bInter)
        %For FO and FOCE WITH interaction, linearize around eta = eta^
        eps = zeros(size(tdata,1),size(sigma,2));
        h = LinMatrixH(model,errormodel,tdata,cdata,theta,eta,eps); %The covariance for this individual
        res_var = diag(diag(h*sigma*h'));
        lC=chol(res_var,'upper')\eye(length(res_var));
        det_res_var = det(res_var);
    else
        %For FO and FOCE WITHOUT interaction, linearize around eta = 0
      %  h = LinMatrixH(tdata,cdata,theta,zeros(size(eta)),eps); %The covariance for this individual
    end

    R=(res'*lC);
    li = -1/2*log(det_res_var)-1/2*(R*R'); % + const
    
else
    %%UDD likelihood
    li=sum(model(tdata,cdata,theta,eta));
end
end