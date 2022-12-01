function [trLabels,tmpMotion2,motion] = mergeClustersByMotion(trLabels,Atr)

labels = unique(trLabels);
motion = zeros(length(labels), length(labels));
countLabel = ones(length(labels),length(labels));
labelNum = length(labels);

% input: labelNum trLabels Atr
% output: motion countLabel
% Atr = full(Atr);
% [motion] = compute_clusters_motion_mex(trLabels, labelNum, Atr);
for i = 1:length(labels)
    for j = i+1:length(labels)
        L = find(trLabels == i);
        R = find(trLabels == j);
        tmpMotion = 0;
        tmpCount = 0;
        for l = 1:length(L)
            idxL = L(l);
            for r = 1:length(R)
                idxR = R(r);
                tmpMotion = tmpMotion + Atr(idxL, idxR);
                if Atr(idxL,idxR)~=0
                    tmpCount = tmpCount + 1;
                end
            end
        end
        motion(i,j) = tmpMotion;%./(length(L)*length(R));
        if tmpCount
            countLabel(i,j) = tmpCount;
        end
    end
end

countLabel1 = ones(labelNum,labelNum);
for i = 1:labelNum
    for j = i+1:labelNum
        inum = length(find(trLabels==i));
        jnum = length(find(trLabels==j));
        countLabel1(i,j) = inum+jnum;
    end
end
% countLabel1 = [countLabel1; zeros(1,labelNum)];

tmpMotion1 = motion./countLabel1;
tmpMotion2 = tmpMotion1>1;
tmpMotion2 = tmpMotion1.*tmpMotion2;
tmpMotion2 = tmpMotion2 + tmpMotion2' + diag(zeros(1,labelNum));

mergeLabels = zeros(1,labelNum);
list = [];
listnum = 1;
list(1) = 1;
l = 1;
lbl = 1;

for i = 1:labelNum
    tmpMotion2(i,:) = [zeros(1,i) tmpMotion2(i,i+1:end)];
end
tmpMotion2 = tmpMotion2 + tmpMotion2';

while(1)
    isfull = 0;
    while ~(l>listnum)
        idx = list(1,l);
        if ~mergeLabels(idx)
            mergeLabels(idx) = lbl;
            % the same row
            tmpL = find(tmpMotion2(idx,:));
            list = [list(1,1:listnum) tmpL];
            listnum = listnum + length(tmpL);
            % the same col
            
        end
        l = l+1;
    end
    isfull = 1;
    for h = 1:labelNum
        if ~mergeLabels(h)
            listnum = listnum+1;
            list(listnum) = h;
            lbl = lbl+1;
            isfull = 0;
            break;
        end
    end
    if(isfull)
        break;
    end
end

for i = 1:labelNum
    tmpidx = find(trLabels == i);
    trLabels(tmpidx) = mergeLabels(1,i);
end

end