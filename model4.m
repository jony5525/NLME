%1-compartment continuous data model at SS
function y=model4(tdata,cdata,theta,eta)

%Model to test case 4, inter occasion variability
%
%tdata = [ID,DV,DV,AGE,RACE,HT,HCTZ,AMT,TIME,EVID,SS,II,OCC,.....], cdata = [COV1,COV2,....,COVn]
tau=12;
dose=1000;

OCC2=1-tdata(:,13); %OCC2 time dependent variable
R=0;
if (tdata(1,5)>=2) %RACE
   R=1;
end

TMP=max((tdata(1,6)-160),0)*theta(1); %HT
TMP2=max(tdata(1,4)-60,0)*theta(2); %AGE
TVCL=theta(3)+TMP-TMP2;
TVCL=TVCL+theta(4)*R+theta(5)*tdata(1,7); %HCTZ
%IF (TVCL.LE.0.) EXIT 1 100
TVV=theta(6)+tdata(1,7)*theta(7);
%IF (TVV.LE.0.) EXIT 1 200
TVKA=theta(8);
cl=TVCL.*exp(eta(4).*tdata(:,13)+eta(7)*OCC2+eta(1));
v =TVV .*exp(eta(5).*tdata(:,13)+eta(8)*OCC2+eta(2));
ka=TVKA.*exp(eta(6).*tdata(:,13)+eta(9)*OCC2+eta(3));
ke=cl./v;

% Steady-State solution of one compartment oral absorption
y=dose./v.*ka./(ka-ke).*(exp(-ke.*tdata(:,9))./(1-exp(-ke*tau)) - exp(-ka.*tdata(:,9))./(1-exp(-ka*tau)));

end