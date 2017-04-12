% derivative of error model w.r.t epsilon
% size of return is (samples per individual x number of sigma)
function h=LinMatrixH(model,errmodel,tdata,cdata,theta,eta,eps)

num_sig = size(eps,2);
%Step length of derivative
h_e = 1E-05;

h = zeros(size(tdata,1),num_sig);

%Complex approximation of derivative of error model w.r.t. epsilon
for k=1:num_sig
    eps_plus=eps;
    eps_plus(:,k) = eps_plus(:,k)+i*h_e;
    err_model_plus = errmodel(model,tdata,cdata,theta,eta,eps_plus);
    h(:,k)=imag(err_model_plus)/h_e;
end
end
