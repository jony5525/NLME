function y=errmodel1(model,tdata,cdata,theta,eta,eps)

%Error model to test case 1
ipred = model(tdata,cdata,theta,eta);

y = ipred.*(1+eps(:,1)); %Proportional residual error

end