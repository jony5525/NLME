% size of return is (samples per individual x (number of sigma x number of omega))
% derivative of model w.r.t. sigma then eta, eval at eps=0 and eta = eta^
function lh=LinMatrixLH(model,errmodel,tdata,cdata,theta,eta,eps)

%Step length of derivative
h_lh = 1E-08;

num_omega =size(eta,1);
num_sigma = size(eps,2);

lh = zeros(size(tdata,1),num_omega*num_sigma);

%Central difference
for k=1:num_omega
    eta_plus=eta;
    eta_minus=eta;
    eta_plus(k) = eta_plus(k)+h_lh;
    eta_minus(k)= eta_minus(k)-h_lh;
    lin_plus  = LinMatrixH(model,errmodel,tdata,cdata,theta,eta_plus,eps);
    lin_minus = LinMatrixH(model,errmodel,tdata,cdata,theta,eta_minus,eps);
    temp = (lin_plus-lin_minus)/(2*h_lh);
    lh(:,(k-1)*num_sigma+1:k*num_sigma)=temp(:,1:num_sigma);
end
end
