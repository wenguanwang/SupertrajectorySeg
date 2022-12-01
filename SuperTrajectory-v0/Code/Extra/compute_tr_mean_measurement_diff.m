function A=compute_tr_mean_measurement_diff(RGB,tr_id,T,sigma,pi,pj,ntr)
tr_start=-1*ones(ntr,1);
tr_end=-1*ones(ntr,1);
[tr_id_on,Is1]=unique(tr_id,'first');
[~,Is2]=unique(tr_id,'last');
tr_start(tr_id_on)=Is1;
tr_end(tr_id_on)=Is2;
pi=uint32(pi);
pj=uint32(pj);
A = compute_tr_mean_measurement_diff_mex(RGB, T,...
    tr_start'-1, tr_end'-1, ntr,...
    sigma, pi-1, pj-1, size(pi,1));
end