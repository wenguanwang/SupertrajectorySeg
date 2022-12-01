function [tr_pos] = get_trajectory_location(tr)
    trNum = length(tr);
    tr_pos = zeros(trNum,3);
    for i = 1:trNum
        each_tr = tr(i);
        pos = each_tr.XYTRGB(1:3,1);
        tr_pos(i,:) = pos';
    end
end