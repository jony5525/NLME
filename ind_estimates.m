%Get the individual empirical bayes estimate for one individual
function eta = ind_estimates(model,errormodel,tdata,cdata,theta,omega,sigma,start_eta,bInter,bUDDLike)

eta_i = start_eta';
c1 = length(start_eta')/2*log(2*pi);
c2 = 1/2*log(det(omega));
c3 = inv(omega);

if (bInter==false && bUDDLike==false) %Calculate only one variance for all samples
    eps = zeros(size(tdata,1),size(sigma,2));
    h = LinMatrixH(model,errormodel,tdata,cdata,theta,zeros(size(eta_i)),eps);
    res_var = diag(diag(h*sigma*h'));
    det_res_var = det(res_var);
    lC=chol(res_var,'upper')\eye(length(res_var)); %Calculate cholesky factorization
else
    lC = 0;
    det_res_var = 0;
end

%eta = fminsearch(@(eta_z) min_function(model,errormodel,tdata,cdata,theta,omega,sigma,bInter,bUDDLike,c1,c2,c3,lC,det_res_var,eta_z),eta_i,optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',50000,'MaxIter',10000));
eta = fminunc(@(eta_z) min_function(model,errormodel,tdata,cdata,theta,omega,sigma,bInter,bUDDLike,c1,c2,c3,lC,det_res_var,eta_z),eta_i,optimset('TolX',1E-10,'TolFun',1E-10,'MaxFunEvals',10000,'Display','Off','LargeScale','off'));
%disp(eta)
end

%This is the function that should be minimized, w.r.t eta
function [ret, dret]=min_function(model,errormodel,tdata,cdata,theta,omega,sigma,bInter,bUDDLike,c1,c2,c3,lC,det_res_var,eta)
    bAnalytic = false;
    li = ind_likelihood(model,errormodel,tdata,cdata,theta,sigma,eta,bInter,lC,det_res_var,bUDDLike); %Individual log likelihood
    ret =c1+c2+1/2*eta*c3*eta'-li;
    if (nargout>1)
        if (~bAnalytic)
            dret=zeros(length(omega),1);
            heta=1e-8;
            for i=1:length(omega)
                eta_plus=eta;
                eta_plus(i)=eta(i)+heta;
                li_plus = ind_likelihood(model,errormodel,tdata,cdata,theta,sigma,eta_plus,bInter,lC,det_res_var,bUDDLike); %Individual log likelihood
                dret(i)=(li_plus-li)/heta;
            end
        else
            dret = sum(dmodel_deta(tdata,cdata,theta,eta));
        end
        dret=c3*eta-dret;
    end
end