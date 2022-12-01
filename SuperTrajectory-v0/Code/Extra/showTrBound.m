function showTrBound(tr, IDX, boundTrLabels, imname, lLen)
    % show the boundary trajectories
    tmpTrColor = zeros(length(tr),1);
    tmpTrLabels = [1:length(tr)];
    tmpLabel1 = boundTrLabels(IDX==1);
    tmpLabel2 = boundTrLabels(IDX==2);
    tmpTrColor(tmpLabel1) = 0.5;
    tmpTrColor(tmpLabel2) = 1;
    h = plot_trajecory_labels_geodist(tr(1:lLen), tmpTrLabels(1:lLen), imname, '_bound', tmpTrColor, 2, [], 1, 4);
end
