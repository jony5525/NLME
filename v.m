function [ret,L] = v(model,errmodel,tdata,cdata,theta,omega,sigma,eta,eps,bInter)

% Covariance of model
% number of samples X number of samples (per individual)

l=LinMatrixL(model,errmodel,tdata,cdata,theta,eta,eps);

if (bInter)
    h=LinMatrixH(model,errmodel,tdata,cdata,theta,eta,eps);
    lh=LinMatrixLH(model,errmodel,tdata,cdata,theta,eta,eps);
    ret = l*omega*l' + diag(diag(h*sigma*h')) + diag(diag(lh*kron(omega,sigma)*lh'));
    L=l;
else
    h=LinMatrixH(model,errmodel,tdata,cdata,theta,zeros(size(eta)),eps);
    ret = l*omega*l' + diag(diag(h*sigma*h'));
    L=l;
end
end
