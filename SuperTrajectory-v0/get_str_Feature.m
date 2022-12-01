function [str_tr, str_co, str_mo, str_lo] = get_str_Feature(tr_lo, tr_mo, tr_co, all_cl)
% trNum = length(tr);
cNum = max(all_cl);
% tr_lens=get_tr_lengths(tr);
% [tr_XYTRGB,tr_id] = quick_tr(tr);
str_lo = zeros(cNum,3);
str_mo = zeros(cNum,2);
str_co = zeros(cNum,3);
str_tr = cell(cNum,1);
for i = 1:cNum
    str_tr{i} = find(all_cl==i);
%     tr_xytrgb = cellfun(@(x){mean(x,2)}, struct2cell(tr(str_tr{i})));
%     tr_xytrgb = cell2mat(tr_xytrgb(1:end));
    co = tr_co(str_tr{i},:);
    str_co(i,:) = mean(co);
    mo = tr_mo(str_tr{i},:);
    str_mo(i,:) = mean(mo);
    lo = tr_lo(str_tr{i},:);
    str_lo(i,:) = mean(lo);
end

