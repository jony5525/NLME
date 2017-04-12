%Individual contributions to the pop likelihood with the FOCE approximation
function ofv_i = ofv_foce(model,errormodel,tdata,cdata,theta,omega,sigma,eta,bInter)

eps = zeros(size(tdata,1),size(sigma,1));

[cov,l] = v(model,errormodel,tdata,cdata,theta,omega,sigma,eta',eps,bInter); %FOCE approximated covariance

tmp_chol = chol(cov,'upper');
lC=tmp_chol\eye(length(cov));
ipred = model(tdata,cdata,theta,eta); %Individual predictions
res = tdata(:,3)-ipred+l*eta';
R=(res'*lC);
ofv_i = log(det(cov))+R*R'; %+ const


end