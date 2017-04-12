function [i_tdata,i_cdata,ind_t]=get_ind_data(start_index,ind,tdata,cdata)
    ind_t = find(tdata(:,1)==ind, 1, 'last');
    i_tdata = tdata(start_index:ind_t,:);
    i_cdata = cdata(ind,:);
end