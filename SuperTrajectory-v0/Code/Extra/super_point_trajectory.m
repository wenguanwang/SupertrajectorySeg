%%如果单独的做点轨迹，应该可以用这些参数来评价点轨迹生成的好坏；因为同一条点轨迹的不同点对应的纹理信
%息应该差别不大。光流是按照块来匹配的，因此用块状超像素来做点超级点轨迹是合理的。在做的时候还可以考虑
%超像素的纹理信息。
function [tr, tr_lo, tr_mo, tr_co, all_center,all_cl]=super_point_trajectory(data, options)
%% Compute Trajectories;
disp('compute trajectories');
% if( exist([options.datafolder '\trajectory.mat'], 'file' ) )
%     load( [options.datafolder '\trajectory.mat'] );
% else  
    [tr] = computeTrLDOF(data, options);
%     save( [options.datafolder '\trajectory.mat'], 'tr', 'trC', '-v7.3' );
% end

[tr_lo, tr_mo, tr_co] = get_trFeature(tr);

%% Compute Trajectory Graph
lens = get_tr_lengths(tr);
t_step = round((data.nframe)/mean(lens));
% t_step = round(mean(lens));
s_step = round(max(sqrt(data.width*data.height/(options.str_num)),2*options.sample_step));
    %away from 2*cuttoff alays affinity zero
options.cutoffx = min(2*s_step,data.width);                  %% spatial neighbour (x)
options.cutoffy = min(2*s_step,data.height);                %% spatial neighbour (y)

% disp('compute trajectory affinities');
% if( exist([options.datafolder '\tra_aff.mat'], 'file' ) )
%     load( [options.datafolder '\tra_aff.mat'] );
% else  
%     [Atr]=computeTrAffinities(tr,trC,options);
%     save( [options.datafolder '\tra_aff.mat'], 'Atr', '-v7.3' );
% end

% inds_keep=sum(Atr,2)>0;
% Atr=Atr(inds_keep,inds_keep);
% tr=tr(inds_keep);
% trC=trC(inds_keep);   

options.verbose = 0;
if options.verbose
    step=3;
    visualize_trajectory_affinity_graph(tr(1:step:end),[],...
        data,Atr(1:step:end,1:step:end),[],[],0,1);
end

%% Compute Super-Trajectory
disp('compute super-trajectory');
[tr_ini_pos] = get_trajectory_location(tr);
% [all_center, all_cl] = get_super_Tra(Atr,tr_ini_pos,data,t_step,s_step,options,tr_lo, tr_ini, tr_end, tr_mo, tr_co, max_v)
[all_center, all_cl] = get_super_Tra(tr,tr_ini_pos,data,t_step,s_step,options,tr_lo,tr_mo, tr_co);

end

