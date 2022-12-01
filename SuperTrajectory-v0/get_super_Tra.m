function [all_center, all_cl] = get_super_Tra(tr,tr_ini_pos,data,t_step,s_step,options,tr_lo,tr_mo, tr_co)
tr_num = size(tr_ini_pos,1);
per_str_num = tr_num/(options.str_num*t_step);
num_thres = 0.2*per_str_num;
all_feature = [tr_lo(:,1:2) tr_co tr_mo];
all_cl = zeros(tr_num,1);
off = 0;
all_center = [];

% all_rows = [];
% all_cols = [];
% for x = 1:s_step:data.width
%     for y = 1:s_step:data.height
%         spixel = find(tr_ini_pos(:,1)>= (x-0.5*s_step)&tr_ini_pos(:,1)<(x+1.5*s_step)&tr_ini_pos(:,2)>=(y-0.5*s_step)&tr_ini_pos(:,2)<(y+1.5*s_step)); 
%         
%         trFrame=get_help_struct(tr(spixel));
%         [pi, pj] = compute_Atraj_neighbour(tr(spixel), trFrame, options.cutoffx, options.cutoffy, options.sample_step,options.verbose);
%         [rows,cols] = pipjtorowcol(pi,pj);
%             
%         all_rows = [all_rows;spixel(rows)];
%         all_cols = [all_cols;spixel(cols)];
%     end
% end
% 
% connet = sparse(all_rows,all_cols,1,tr_num,tr_num);

for x = 1:s_step:data.width
    for y = 1:s_step:data.height
        spixel = find(tr_ini_pos(:,1)>=x&tr_ini_pos(:,1)<(x+s_step)&tr_ini_pos(:,2)>=y&tr_ini_pos(:,2)<(y+s_step));
%         if size(spixel,1)>=num_thres
%             dist = Atr(spixel,spixel);
%             connect = sum(dist);
%             spixel(connect==0) = [];
%         end
        if size(spixel,1)<num_thres
            continue;
        else
            trFrame=get_help_struct(tr(spixel));
            [pi, pj] = compute_Atraj_neighbour(tr(spixel), trFrame, options.cutoffx, options.cutoffy, options.sample_step,options.verbose);
            [rows,cols] = pipjtorowcol(pi,pj);
            
            dist2 = exp(-(EuclideanDistance(tr_lo(spixel,1:2),tr_lo(spixel,1:2)) + ...
            EuclideanDistance(tr_co(spixel,:),tr_co(spixel,:)) + ...
        EuclideanDistance(tr_mo(spixel,:),tr_mo(spixel,:))));
            dist = zeros(size(dist2));
            inx = sub2ind(size(dist2),rows(:),cols(:));
            dist(inx) = dist2(inx);
            dist(logical(speye(size(dist)))) = 0;
%             dist = Atr(spixel,spixel);
            center_num = ceil(size(spixel,1)/per_str_num);
            [center_ins,cl]= cluster_dp(dist,center_num);
            all_center = [all_center;spixel(center_ins)];
            all_cl(spixel)=cl+off;
            off = off+max(cl(:));
        end
    end
end

ins_off = find(all_cl==0);

% dist = Atr(ins_off,all_center);

% query_feature = all_feature(ins_off,:);
% test_feature = all_feature(all_center,:);
% dist = exp(-(EuclideanDistance(query_feature,test_feature) ));
if ~isempty(ins_off)
    all_ins = [ins_off;all_center];
    trFrame=get_help_struct(tr([ins_off;all_center]));
    [pi, pj] = compute_Atraj_neighbour(tr([ins_off;all_center]), trFrame, options.cutoffx, options.cutoffy, options.sample_step,options.verbose);
    [rows,cols] = pipjtorowcol(pi,pj);

    dist2 = exp(-(EuclideanDistance(tr_lo(ins_off,1:2),tr_lo(all_center,1:2)) + ...
                EuclideanDistance(tr_co(ins_off,:),tr_co(all_center,:)) + ...
            EuclideanDistance(tr_mo(ins_off,:),tr_mo(all_center,:))));
        
    dist = zeros(size(all_ins,1),size(all_ins,1));
    inx = sub2ind(size(dist),rows(:),cols(:));
    dist(inx) = 1;
    dist(logical(speye(size(dist)))) = 0;
    dist = dist(1:size(ins_off,1),end-size(all_center,1)+1:end);
    dist = dist.*dist2;

    [max_d,max_ins] = max(dist,[],2);
    all_cl(ins_off(max_d>0)) = max_ins(max_d>0);
    all_cl(ins_off(max_d==0)) = max(all_cl(:))+[1:size(ins_off(max_d==0),1)];
    all_center = [all_center;ins_off(max_d==0)];
end
% dist = Atr;
% dist = dist(:,all_center);
% [max_d,max_ins] = max(dist,[],2);

query_feature = all_feature;
test_feature = all_feature(all_center,:);
[sort_Ind,~] = KDSearch(single(test_feature'),single(query_feature'),1,1);
max_ins = sort_Ind';
    
all_cl2 = max_ins;
all_cl2(all_center) = all_cl(all_center);
all_cl = all_cl2;

for loop=1:10
  loop  
all_cl_new = zeros(tr_num,1);
off = 0;
all_center_new = [];
for i = 1:size(all_center,1)
     spixel = find(all_cl == i);
     if size(spixel,1)<num_thres
        continue;
     else 
%         dist = Atr(spixel,spixel);
        trFrame=get_help_struct(tr(spixel));
        [pi, pj] = compute_Atraj_neighbour(tr(spixel), trFrame, options.cutoffx, options.cutoffy, options.sample_step,options.verbose);
        [rows,cols] = pipjtorowcol(pi,pj);

        dist2 = exp(-(EuclideanDistance(tr_lo(spixel,:),tr_lo(spixel,:)) + ...
            EuclideanDistance(tr_co(spixel,:),tr_co(spixel,:)) + ...
        EuclideanDistance(tr_mo(spixel,:),tr_mo(spixel,:))));
    
        dist = zeros(size(dist2));
        inx = sub2ind(size(dist2),rows(:),cols(:));
        dist(inx) = dist2(inx);
        dist(logical(speye(size(dist)))) = 0;
       
        center_num = 1;
        [center_ins,cl]= cluster_dp(dist,center_num);
        all_center_new = [all_center_new;spixel(center_ins)];
        all_cl_new(spixel)=cl+off;
        off = off+max(cl(:));
     end
end


ins_off = find(all_cl_new==0);
if ~isempty(ins_off)
    all_ins = [ins_off;all_center];
    trFrame=get_help_struct(tr([ins_off;all_center]));
    [pi, pj] = compute_Atraj_neighbour(tr([ins_off;all_center]), trFrame, options.cutoffx, options.cutoffy, options.sample_step,options.verbose);
    [rows,cols] = pipjtorowcol(pi,pj);

    dist2 = exp(-(EuclideanDistance(tr_lo(ins_off,1:2),tr_lo(all_center,1:2)) + ...
                EuclideanDistance(tr_co(ins_off,:),tr_co(all_center,:)) + ...
            EuclideanDistance(tr_mo(ins_off,:),tr_mo(all_center,:))));
        
    dist = zeros(size(all_ins,1),size(all_ins,1));
    inx = sub2ind(size(dist),rows(:),cols(:));
    dist(inx) = 1;
    dist(logical(speye(size(dist)))) = 0;
    dist = dist(1:size(ins_off,1),end-size(all_center,1)+1:end);
    dist = dist.*dist2;

    [max_d,max_ins_new] = max(dist,[],2);
    all_cl_new(ins_off(max_d>0)) = max_ins_new(max_d>0);
    all_cl_new(ins_off(max_d==0)) = max(all_cl_new(:))+[1:size(ins_off(max_d==0),1)];
    all_center_new = [all_center_new;ins_off(max_d==0)];
end

query_feature = all_feature;
test_feature = all_feature(all_center_new,:);
[sort_Ind,~] = KDSearch(single(test_feature'),single(query_feature'),1,1);
max_ins = sort_Ind';

all_cl2 = max_ins;
all_cl2(all_center_new) = all_cl_new(all_center_new);
all_cl_new = all_cl2;

all_center = all_center_new;
all_cl = all_cl_new; 
end




end