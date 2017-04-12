%Individual contributions to the pop likelihood with the FOCEI nonmem approximation
function ofv_i = ofv_focei_nonmem(model,errormodel,tdata,cdata,theta,omega,sigma,eta)

eps = zeros(1,size(sigma,1));

l=LinMatrixL(model,errormodel,tdata,cdata,theta,eta',eps);

h=LinMatrixH(model,errormodel,tdata,cdata,theta,eta',eps);

cov1 = diag(diag(h*sigma*h'));

inv_cov1 = inv(cov1);

inv_om = inv(omega);

li = ind_likelihood(model,errormodel,tdata,cdata,theta,sigma,eta,true,0,0,false);
hess=0;
for i=1:size(tdata,1)
  tmp=deta_dv(model,errormodel,tdata(i,:),cdata,theta,sigma,eta,eps); 
  hess=hess+l(i,:)'*inv_cov1(i,i)*l(i,:)+1/2*inv_cov1(i,i)*tmp'*inv_cov1(i,i)*tmp;
end
       
ofv_i = -2*li+log(det(omega))+eta*inv_om*eta'+log(det(inv_om+hess));


end

function grad = deta_dv(model,errormodel,tdata,cdata,theta,sigma,eta,eps)
grad=zeros(size(eta));
h=1E-05;
for i=1:length(eta)
    eta_plus = eta;
    eta_minus =eta;
    eta_plus(i) = eta_plus(i) + h;
    eta_minus(i) = eta_minus(i) - h;
    hplus=LinMatrixH(model,errormodel,tdata,cdata,theta,eta_plus,eps);
    hminus=LinMatrixH(model,errormodel,tdata,cdata,theta,eta_minus,eps);
    ff_plus = diag(diag(hplus*sigma*hplus'));
    ff_minus= diag(diag(hminus*sigma*hminus'));
    grad(i) = (ff_plus-ff_minus)/(2*h);
end

end


