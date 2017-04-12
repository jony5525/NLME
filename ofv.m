function [ofv_sum] = ofv(model,errormodel,type,tdata,cdata,theta,omega,sigma,ebes,bInter,bUDDLike,bReport)

if (bUDDLike && type<3)
    error('User defined ln likelihood is only available with the Laplace or MC methods!');
end

num_ind = max(tdata(:,1));
start_index=1;

n_omega = size(omega,1);
ofv_sum = 0;


for i=1:num_ind
    [i_tdata,i_cdata,end_index] = get_ind_data(start_index,i,tdata,cdata);
    start_index = end_index+1;
    eta = zeros(n_omega,1);
    switch (type)
        case 0 %First order approximation (FO)
            ofv_i = ofv_fo(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,bInter);
        case 1 %First order conditional approximation (FOCE)
            eta = ind_estimates(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,bInter,bUDDLike);
            ofv_i = ofv_foce(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,bInter);
        case 2 % FOCEI (NONMEM Way)
            eta = ind_estimates(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,true,bUDDLike);
            ofv_i = ofv_focei_nonmem(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta);
        case 3 %Laplace approximation
            eta = ind_estimates(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,bInter,bUDDLike);
            ofv_i=ofv_laplace(model,errormodel,i_tdata,i_cdata,theta,omega,sigma,eta,bInter,bUDDLike);
    end
    
    ofv_sum = ofv_sum+ofv_i;
    if (bReport)
        fprintf('The -2ll for individual %i is: %d\n',i,ofv_i);
    end
end

end
