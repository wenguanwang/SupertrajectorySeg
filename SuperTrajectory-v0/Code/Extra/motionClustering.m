function [tr_labels_rot,geoNormDist,foreSeeds,fore,foreNormDist]=motionClustering(tr,Vtr,Str,imnames,...
    para,verbose,lLen,Atr)
K=size(Vtr,2);
binsoltr = getbinsol(Vtr);
tr_labels_rot = full(sum(binsoltr.*repmat(1:K, size(binsoltr,1),1),2)');
tiny_cls = get_tiny_cls(tr_labels_rot, para);
tr_labels_rot(ismember(tr_labels_rot,tiny_cls))=0;
[tr_labels_rot] = make_label_continuous(tr_labels_rot);
tr_labels_rot(tr_labels_rot<0)=0;
  

[geoNormDist foreNormDist] = mergingClustersByMotion(tr,tr_labels_rot,Atr,imnames,para,lLen);


[foreSeeds fore] = findForeSeeds(tr, geoNormDist, foreNormDist, para, imnames, lLen);
% if ~isempty(foreNormDist)
%     [foreSeedsRes foreRes] = findForeSeeds(tr, foreNormDist, para, imnames, lLen);
% end


% h=plot_trajectory_labels(tr, tr_labels, imnames(5), 2, [], 1, 4);
% h_l=plot_trajectory_labels(tr(1:lLen), tr_labels(1:lLen), imnames(5), '_l', 2, [], 1, 4);
% h_r=plot_trajectory_labels(tr(lLen+1:end), tr_labels(lLen+1:end),imnames(10), '_r', 2, [], 1, 4);


end
