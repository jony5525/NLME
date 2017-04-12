%Individual contributions to the pop likelihood with the second order
%approximation (Laplace in NONMEM)
function [ofv_i] = ofv_laplace(model,errmodel,tdata,cdata,theta,omega,sigma,eta,bInter,bUDDLike)

bAnalytic = false;


if (~bInter && ~bUDDLike)    
    eps = zeros(size(tdata,1),size(sigma,2));
    h = LinMatrixH(model,errmodel,tdata,cdata,theta,zeros(size(eta_i)),eps);
    res_var = diag(diag(h*sigma*h'));
    lC=chol(res_var,'upper')\eye(length(res_var)); %Calculate cholesky factorization
    det_res_var = det(res_var);
else
    lC = 0;
    det_res_var = 0;
end


ofv_i = ind_likelihood(model,errmodel,tdata,cdata,theta,sigma,eta,bInter,lC,det_res_var,bUDDLike);

inv_om = inv(omega);

if (~bAnalytic)
    hess=hessian_eta(model,errmodel,tdata,cdata,theta,sigma,bInter,lC,det_res_var,bUDDLike,eta);
else
    hess = sum(d_dmodel_deta(tdata,cdata,theta,eta));
end

X  = inv_om-hess;
ofv_i = -2*ofv_i + log(det(omega))+eta*inv_om*eta'+log(det(X));
end

function hess=hessian_eta(model,errmodel,tdata,cdata,theta,sigma,bInter,lC,det_res_var,bUDDLike,eta)
n=length(eta); 
hess=zeros(n);   % Memory for the Hessian matrix
h=1E-03;                      
h2=h*h;  %Inner step length  (at least as small as h)                         
for k=1:n                           
    eta_plus = eta;
    eta_minus = eta;

    eta_plus(k) = eta_plus(k)+h;
    eta_minus(k) = eta_minus(k)-h;
    
    g_plus = grad_eta(model,errmodel,tdata,cdata,theta,sigma,bInter,lC,det_res_var,bUDDLike,eta_plus,h2,k);
    g_minus = grad_eta(model,errmodel,tdata,cdata,theta,sigma,bInter,lC,det_res_var,bUDDLike,eta_minus,h2,k);

    g = (g_plus-g_minus)./(2*h);
    hess(k,k:n) = g(k:n);
end
hess = diag(diag(hess))+triu(hess,1) + triu(hess,1)';
end


function grad=grad_eta(model,errmodel,tdata,cdata,theta,sigma,bInter,lC,det_res_var,bUDDLike,eta,h,k)
grad=zeros(size(eta));
for i=k:length(eta)
    eta_plus = eta;
    eta_minus =eta;
    eta_plus(i) = eta_plus(i) + h;
    eta_minus(i) = eta_minus(i) - h;
    ff_plus =ind_likelihood(model,errmodel,tdata,cdata,theta,sigma,eta_plus,bInter,lC,det_res_var,bUDDLike);
    ff_minus=ind_likelihood(model,errmodel,tdata,cdata,theta,sigma,eta_minus,bInter,lC,det_res_var,bUDDLike);
    grad(i) = (ff_plus-ff_minus)/(2*h);
end
end