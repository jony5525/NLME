function l=LinMatrixL(model,errmodel,tdata,cdata,theta,eta,eps)
% Derivative of individul predictions (model without residual error) w.r.t eta
% size: (samples per individual x number of random effects)
%
%Step length of derivative
h_l = 1E-05;

num_omega =size(eta,1);
l=zeros(size(tdata,1),num_omega);

%Complex approximation of derivative of model w.r.t. eta
for k=1:num_omega
    eta_plus=eta;
    eta_plus(k) = eta_plus(k)+i*h_l;
    err_model_plus = errmodel(model,tdata,cdata,theta,eta_plus,eps);
    l(:,k)=imag(err_model_plus)/h_l;
end
end


