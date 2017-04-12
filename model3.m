function y=model3(tdata,cdata,theta,eta)

%Model to test case 3
% 
%tdata = [ID,DV,DOS,.....], cdata = [COV1,COV2,....,COVn]
 
bas= theta(1)*exp(eta(1));
d50 = theta(2);
lam=bas*(1-tdata(:,3)./(tdata(:,3)+d50));

lnfac = (tdata(:,2)>1).*(tdata(:,2).*log(tdata(:,2)+1E-16)-tdata(:,2)+log((tdata(:,2)+1E-16).*(1+4*tdata(:,2).*(1+2*tdata(:,2))))/6+log(3.1415)/2);
 
y=-lam+tdata(:,2).*log(lam)-lnfac;

end