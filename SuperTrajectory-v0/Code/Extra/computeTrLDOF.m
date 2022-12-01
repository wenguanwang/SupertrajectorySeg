function [tr] = computeTrLDOF(data, options)
%% link flow fields
[tr] = linkFlowLDOF(data, options);
lens=get_tr_lengths(tr);
tr=tr(lens > 3);
% trC=trC(lens > 5);
end