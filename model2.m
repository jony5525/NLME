function y=model2(tdata,cdata,theta,eta)

%Model to test case 2

%tdata = [ID,TIME,DV,.....], cdata = [COV1,COV2,....,COVn]

V = theta(1)*exp(eta(1));
KE = theta(2)*exp(eta(2));
KA = theta(3)*exp(eta(3));

y = 320/V*KA/(KA-KE).*(exp(-KE.*tdata(:,2))-exp(-KA.*tdata(:,2)));

end