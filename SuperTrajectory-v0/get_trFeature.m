function [tr_lo, tr_mo, tr_co] = get_trFeature(tr)
% trNum = length(tr);
% [tr_XYT,tr_id] = quick_tr(tr);
% max_v = zeros(data.nframe,1);
% tr_co = zeros(trNum,3);
% tr_ini = zeros(trNum,1);
% tr_end = zeros(trNum,1);
% tr_lo = zeros(trNum,3);
% tr_mo = zeros(trNum,2);
aggr = 3;
tr_xytrgb = cellfun(@(x){mean(x,2)}, struct2cell(tr));
tr_xytrgb = cell2mat(tr_xytrgb(1:end));
tr_lo = tr_xytrgb(1:3,:)';
tr_co = tr_xytrgb(4:6,:)';
tr_mo = cellfun(@(x){mean((x(1:2,aggr:end)-x(1:2,1:end-aggr+1)),2)}, struct2cell(tr));
tr_mo = cell2mat(tr_mo(1:end))';
% for i = 1:trNum
%     tr_ALLlo = tr_XYT(:,tr_id==i)';
%     tr_loind = sub2ind([data.height data.width data.nframe], int16(tr_ALLlo(:,2)), int16(tr_ALLlo(:,1)), int16(tr_ALLlo(:,3)));
%     
%     tr_R = data.AllframesR(tr_loind);
%     tr_G = data.AllframesG(tr_loind);
%     tr_B = data.AllframesB(tr_loind);    
%     tr_co(i,:) = mean(double([tr_R tr_G tr_B]));
%     
%     tr_ini(i) = tr_ALLlo(1,3);
%     tr_end(i) = tr_ALLlo(end,3);
% 
%     tr_lo(i,:) = mean(tr_ALLlo);
%     
%     v=tr_ALLlo(aggr:end,1:2)-tr_ALLlo(1:end-aggr+1,1:2);
%     x = sum(v.^2,2)+1;
%     max_v(tr_ALLlo(1:end-aggr+1,3)) = max([max_v(tr_ALLlo(1:end-aggr+1,3)) x],[],2);
%     
%     tr_mo(i,:) = mean(v);
% end
% max_v(max_v<1) = 1;
% max_v = sqrt(max_v);

 
