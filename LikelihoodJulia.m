%%%%% Evaluation of OBJ (-2 log likelihood) for different approximation
%%%%% methods for NLMEM (Non Linear Mixed Effect Models)

%%%%% Examples based on code from Joakim Nyberg joakim.nyberg@farmbio.uu.se
%%%%% The different OFV apporximations are written by Joakim Nyberg but some of
%%%%% the code;
%%%%% (EBE estimation, LinMatrixL, LinMatrixH, V) comes from the PopED sofware written by Andrew
%%%%% Hooker (andrew.hooker@farmbio.uu.se),
%%%%% Joakim Nyberg (joakim.nyberg@farmbio.uu.se) and
%%%%% Sebastian Ueckert (sebastian.ueckert@farmbio.uu.se)

%%%%% Likelihoods are assuming that parameters are described as point
%%%%% estimates from distributions (e.g. similar to the NONMEM software)
%%%%% as opposed to distributions with link functions (e.g. similar to the
%%%%% MONOLIX software).

%%%%% References:
%%%%% M Davidian, D Giltinan: Nonlinear Models for Repeated Measurements data,
%%%%% Chapman and Hall/CRC, 1989
%%%%% Wang, Y. J. Derivations of various NONMEM estimation methods
%%%%% Pharmacokinet Pharmacodyn (2007) 34: 575. doi:10.1007/s10928-007-9060-6
%%%%% M Lavielle, Mixed Effects Models for the Population Approach: Models, Tasks, Methods and Tools
%%%%% July 14, 2014 by Chapman and Hall/CRC


%%%%% First example evaluation simple continuous data model with FO, FOI, FOCE, FOCEI and LAPI
%%%%% Three parameters; one typical value, THETA, one
%%%%% Inter-Individual-Variance (IIV) OMEGA and one residual variance
%%%%% (SIGMA)

%Estimation type 0=FO, 1 = FOCE, 2 = FOCEI NONMEM WAY, 3 = Laplace, 4 = MC
etype = 0;
%Interaction or not
bInter = 0;
%UDD likelihood or "normal" likelihood
bUDDLike = false;

%Report = true (print), false (silent)
bReport = true;

%Simulate similar data as the model (same design, model and parameters)
%(true or false)
bSimpleSim = true;

%Fixed effects
theta = 0.5;

%Random IIV effect
omega = 0.04;

%Residual random effect
sigma = 0.1;

%cdata = constant dependent datafile for all individuals
cdata = ones(20,1).*10;
%Read in time-dependent data for all individuals, 
tdata = csvread('sim_data_model1.csv');

errmodel=@errmodel1;
model=@model1;

% ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,1),bInter,bUDDLike,bReport);
% fprintf('The -2ll for model 1 with FO is: %3.15f\n',ofv_sum);
% ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,1),true,bUDDLike,bReport);
% fprintf('The -2ll for model 1 with FO with interaction is: %3.15f\n',ofv_sum);
% etype=1;
% ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,1),bInter,bUDDLike,bReport);
% fprintf('The -2ll for model 1 with FOCE is: %3.15f\n',ofv_sum);
% etype=2;
% ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,1),true,bUDDLike,bReport);
% fprintf('The -2ll for model 1 with FOCE with interaction is: %3.15f\n',ofv_sum);
% etype=3;
% ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,1),true,bUDDLike,bReport);
% fprintf('The -2ll for model 1 with Laplace with interaction is: %3.15f\n',ofv_sum);

%%%%% NONMEM OUTPUT (-2LL)
%%%%% FO:    56.474912258258158
%%%%% FOI:   NA
%%%%% FOCE:  56.476216665029462 
%%%%% FOCEI: 56.410938825140313
%%%%% LAPLI: 56.810343602063618


%%%%% Second example evaluation one compartment continuous data model with FO, FOI, FOCE, FOCEI and LAPI
%%%%% 11 parameters; 3 typical values, THETAs, 3
%%%%% Inter-Individual-Variance (IIV) OMEGA with 3 covariance and two
%%%%% residual variances (SIGMA)

%Fixed effects
theta = [3.24467E+01  8.72879E-02  1.49072E+00];
  
%Random IIV effect
omega = [1.93973E-02 1.20854E-02 5.69131E-02
    1.20854E-02  2.02375E-02 -6.47803E-03
    5.69131E-02 -6.47803E-03  4.34671E-01];

%Residual random effect
sigma = [1.70385E-02 0; 0 8.28498E-02];

%cdata = constant dependent datafile for all individuals
cdata = ones(20,1).*12;
%Read in time-dependent data for all individuals, 
tdata = csvread('THEOPP.csv');

errmodel=@errmodel2;
model=@model2;

etype=0;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),false,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FO is: %3.15f\n',ofv_sum);
etype=0;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FO with interaction is: %3.15f\n',ofv_sum);
etype=1;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),false,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FOCE is: %3.15f\n',ofv_sum);

etype=2;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FOCE with interaction is: %3.15f\n',ofv_sum);
etype=3;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 2 with LAPLACE with interaction is: %3.15f\n',ofv_sum);

%%%%% NONMEM OUTPUT (-2LL)
%%%%% FO:    105.45097989219278
%%%%% FOI:   NA
%%%%% FOCE:  94.838182670072669
%%%%% FOCEI: 92.830868885461669
%%%%% LAPLI: 93.641663857694311

%%%%% Third example evaluation of count (Poisson) data with LAPLI
%%%%% 3 parameters; 2 typical values, THETAs, 1
%%%%% Inter-Individual-Variance (IIV) OMEGA
%%%%% A user defined individual likelihood for Possion data is specified as
%%%%% a log likelihood in the model

%User defined individual likelihood
bUDDLike=true;

%Fixed effects
theta = [1.02930E+00  4.51854E-01];

%Random IIV effect
omega = 1.22009E-01;

%Residual random effect
sigma = [];

%tdata = time dependent datafile for all individuals, id,dv,dose,....
tdata = csvread('sim_poisson.csv');
%Individual constant data (not used)
cdata = zeros(20,1);

model=@model3;

%Note - No residual error model for non-continuous data

etype=3;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omega,sigma,size(omega,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 3 with LAPLACE with interaction is: %3.15f\n',ofv_sum);

%%%%% NONMEM OUTPUT (-2LL)
%%%%% LAPLI: 3809.8059929876335

%%%%% Fourth example evaluation one compartment continuous data model with FO, FOI, FOCE, FOCEI, LAPL and LAPI and inter occasion variability
%%%%% 11 parameters; 8 typical values, THETAs, 8
%%%%% Inter-Individual-Variance (IIV) OMEGA with 3 covariance and two
%%%%% Inter-Occasion-Variance (IOV) OMEGA with 3 covariance and
%%%%% residual variances (SIGMA)

%User defined individual likelihood
bUDDLike=false;

%Fixed effects
theta = [4.00000E-01 3.00000E-01 2.30000E+01 1.90000E+00 -4.80000E+00 1.03000E+02 -2.70000E+01 1.40000E+00];
  
%Random IIV effect
omega = [0.2 0.1 0
         0.1 0.2 0
         0   0   1.2];

%Random IOV effect
gamma = [0.2 0.1 0
         0.1 0.2 0
         0   0   0.5];

omeganew = blkdiag(omega,gamma,gamma); %Make a big omega-matrix with 2 occasions    
     

%Residual random effect
sigma = [0.1];

%cdata = constant dependent datafile for all individuals
cdata = ones(64,1).*0;
%Read in time-dependent data for all individuals,
tdata = csvread('praz21_sim_matlab1.csv');

errmodel=@errmodel4;
model=@model4;

etype=0;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omeganew,2),false,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FO is: %3.15f\n',ofv_sum);

etype=0;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omeganew,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 2 with FO with interaction is: %3.15f\n',ofv_sum);

etype=1;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omega,2),false,bUDDLike,bReport);
fprintf('The -2ll for model 4 with FOCE is: %3.15f\n',ofv_sum);
% 
etype=2;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omeganew,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 4 with FOCE with interaction is: %3.15f\n',ofv_sum);

etype=3;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omeganew,2),false,bUDDLike,bReport);
fprintf('The -2ll for model 4 with LAPLACE WITHOUT interaction is: %3.15f\n',ofv_sum);

etype=3;
ofv_sum= ofv(model,errmodel,etype,tdata,cdata,theta,omeganew,sigma,size(omeganew,2),true,bUDDLike,bReport);
fprintf('The -2ll for model 4 with LAPLACE with interaction is: %3.15f\n',ofv_sum);

%%%%% NONMEM OUTPUT (-2LL)
%%%%% FO:    3934.9931582677168
%%%%% FOI:   NA
%%%%% FOCE:  3074.1192663684797
%%%%% FOCEI: 949.85682322100604
%%%%% LAPL:  3120.3839170232632
%%%%% LAPLI: 969.53263286478852 (MCETA=100), ~1018.0647 without MCETA


