function [A]=computeTrAffinities(tr, trC, options)
lens = get_tr_lengths(tr);
cutoffx = options.cutoffx;                  %% spatial neighbour (x)
cutoffy = options.cutoffx;                %% spatial neighbour (y)
max_tr_aggr = min(5,lens);
% my_var_euk = options.my_var_euk;
my_var_euk = max(cutoffx,cutoffy);
my_var_cuk = options.my_var_cuk;
my_var_ut = options.my_var_ut;
sample_step = options.sample_step;
verbose = options.verbose;
ntr=length(tr);

%% neighbors
trFrame=get_help_struct(tr);
[pi, pj] = compute_Atraj_neighbour(tr, trFrame, cutoffx, cutoffy, sample_step,verbose);
%tr     :  structure of trajectories. Each index has field XYTPos
%trFrame:  cell array of space time point.
%           Each cell is #Npts x 3.
%                #Npts is number of space time point,
%                 columns are (x, y, frameid)
% cutoffx, cutoffy : where to cutout
% samplestep : how much to sample
% verbose : if >=2, display info
% [pi, pj]        : index pair representation for MATLAB sparse matrix   
%前【pj(n):pj(n+1)】对应着和点轨迹n相连的点轨迹pi【pj(n):pj(n+1)】
[rows,cols] = pipjtorowcol(pi,pj);
clear pi pj

%% spatial distance
[XYT,tr_id]=quick_tr(tr);
RGB = cat( 2, trC.RGB );
fprintf('Compute location distance \n');
vals_euc=compute_tr_max_measurement_diff(XYT(1:2,:)',tr_id,XYT(3,:)',my_var_euk,rows,cols,ntr);
%% color distance
fprintf('Compute color distance \n');
vals_cuc=compute_tr_mean_measurement_diff(RGB',tr_id,XYT(3,:)',my_var_cuk,rows,cols,ntr);

max_aggr=5;
min_aggr=3;

%% velocity
fprintf('Compute motion distance \n');
cnt=0;
% my_var_utc=my_var_ut;
for aggr=max_aggr:-1:min_aggr
    %
    max_v = zeros(size(trFrame,2),1);
    cnt=cnt+1;
    tr_id_on=find(max_tr_aggr>=aggr);
    keep=find(ismember(rows,tr_id_on)& ismember(cols,tr_id_on));
%     fprintf('len:%d',sum(keep));
    if sum(keep)==0
        break;
    end
%     fprintf('%d on trajectories',length(tr_id_on));
    clear Us T tr_ids
    for ii=tr_id_on'
        XYTc=tr(ii).XYTPos;
        v=XYTc(1:2,aggr:end)-XYTc(1:2,1:end-aggr+1);
        v= [v repmat(v(:,end),1,aggr-1)];
        if isempty(v)
            error('error!');
        end
        x = sum(v.^2)+1;
        max_v(XYTc(3,1:end)) = max([max_v(XYTc(3,1:end)) x'],[],2);
        Us{ii}=v;
        T{ii}=XYTc(3,1:end);
        tr_ids{ii}=ii*ones(size(T{ii}));
    end
    max_v(max_v<1) = 1;
    max_v = sqrt(max_v);
    tr_id=cat(2,tr_ids{:})';
    Tall=cat(2,T{:})';
    Usall=cat(2,Us{:})';
    vals=compute_tr_max_measurement_diff_ut(Usall,tr_id,Tall,max_v,rows(keep),cols(keep),ntr);   
%     vals=compute_tr_max_measurement_diff_ut(Usall,tr_id,Tall,my_var_utc,rows(keep),cols(keep),ntr);
%     keep_on=vals>0;
%     Cols{cnt}=cols(keep(keep_on));
%     Rows{cnt}=rows(keep(keep_on));
%     Vals{cnt}=vals(keep_on).*vals_euc(keep(keep_on));
    
    Cols{cnt}=cols(keep);
    Rows{cnt}=rows(keep);
    Vals{cnt}=exp(-(options.m1*vals+options.m2*vals_cuc(keep)+vals_euc(keep)));
    rows(keep)=[];
    cols(keep)=[];
%     my_var_utc=(max_aggr^2/aggr^2)*my_var_ut;
end
A=sparse(cat(1,Rows{:}),cat(1,Cols{:}),cat(1,Vals{:}),ntr,ntr);
end

