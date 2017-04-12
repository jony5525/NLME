function y=model(tdata,cdata,theta,eta)

%Model to test case 1

%tdata = [ID,TIME,DV,.....], cdata = [COV1,COV2,....,COVn]

ke = theta(1)*exp(eta(1));
y = cdata(1)*exp(-ke*tdata(:,2));

end