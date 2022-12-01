function [seg]=seg_str(data, options,tr,all_cl,str_tr, str_lo, str_mo, str_co)
t = 1;
[tr_XYTRGB,tr_id]=quick_tr(tr);
tr_ins = find(tr_XYTRGB(3,:)==t);

str_ins = all_cl(tr_id(tr_ins));
tr_lo = tr_XYTRGB(1:2,tr_ins)';
tr_loind = sub2ind([data.height data.width], int16(tr_lo(:,2)), int16(tr_lo(:,1)));
label = data.gt(tr_loind);
fore_tr = tr_id(tr_ins);
fore_tr(label==0)=[];
back_tr = tr_id(tr_ins);
back_tr(label==1)=[];

num_str = size(str_tr,1);
str_label = -1*ones(num_str,1);

fore_str = str_ins(label);
back_str = str_ins(~label);

str_label(fore_str) = 1;
str_label(back_str) = 1;
un_str = find(str_label==-1);
la_str = find(str_label==1);
% str_W = zeros(num_str,num_str);
% str_W(logical(speye(size(str_W)))) = 1;

str_W = exp(-(...%EuclideanDistance(str_lo(un_str,1:2),str_lo(la_str,1:2)) + ...
            EuclideanDistance(str_co(un_str,:),str_co(la_str,:)) + ...
        EuclideanDistance(str_mo(un_str,:),str_mo(la_str,:))));    
% [sort_dis sort_ins]=sort(str_W ,2,'descend');
% sort_dis = sort_dis(:,1:25);
% sort_ins = sort_ins(:,1:25);
% query_ins = repmat([1:size(un_str,1)]',[1 25]);
% str_W(:) = 0; 
% str_W(sub2ind(size(str_W), query_ins(:), sort_ins(:) )) = sort_dis(:);

Hist_fore = hist(fore_str(:),1:1:num_str);
Hist_back = hist(back_str(:),1:1:num_str);
str_pro = Hist_fore./(Hist_fore+Hist_back+0.001);
str_pro = str_pro(la_str);
% E = sparse(1:num_str,1:num_str,ones(num_str,1)); 
% iD = sparse(1:num_str,1:num_str,1./sum(str_W));
iD = repmat(1./sum(str_W,2),[1 size(la_str,1)]);
% P = iD*str_W;
P = iD.*str_W;
% v = str_pro';
% for T = 1:10
%    v = P*v ;
%    v(str_label==1) = str_pro(str_label==1); 
% end
v = zeros(num_str,1);
v(un_str) = P*(str_pro');
v(la_str) = str_pro;

options.verbose = 0;
if options.verbose
    for t=1:data.nframe
        [tr_XYT,tr_id]=quick_tr(tr);
        tr_ins = find(tr_XYT(3,:)==t);
        str_ins = all_cl(tr_id(tr_ins));
        tr_lo = tr_XYT(1:2,tr_ins)';
        imshow(ones(data.height,data.width,3)*100);
        hold on;
        scatter(int16(tr_lo(:,1)),int16(tr_lo(:,2)),3*ones(size(tr_lo,1),1),v(str_ins,:),'filled');
    end
end

tr_loind = sub2ind([data.height data.width data.nframe], int16(tr_XYTRGB(2,:)), int16(tr_XYTRGB(1,:)), int16(tr_XYTRGB(3,:)));
video_l = zeros(data.height,data.width,data.nframe);
video_label = zeros(data.height,data.width,data.nframe);
% video_l(tr_loind) = str_label(all_cl(tr_id));
video_l(:,:,1) = double(data.gt);
video_label(:,:,1) = 1;
video_l(tr_loind) = v(all_cl(tr_id));
video_label(tr_loind) = 1;

tr_xytrgb = struct2cell(tr(fore_tr));
tr_xytrgb = cell2mat(tr_xytrgb(1:end));
fore_tr_loind = sub2ind([data.height data.width data.nframe], int16(tr_xytrgb(2,:)), int16(tr_xytrgb(1,:)), int16(tr_xytrgb(3,:)));
video_l(fore_tr_loind) = 1;

tr_xytrgb = struct2cell(tr(back_tr));
tr_xytrgb = cell2mat(tr_xytrgb(1:end));
back_tr_loind = sub2ind([data.height data.width data.nframe], int16(tr_xytrgb(2,:)), int16(tr_xytrgb(1,:)), int16(tr_xytrgb(3,:)));
video_l(back_tr_loind) = 0;

[s_fpro] = getSuperpixelFprob( double(video_l), video_label ,data.superpixels, data.slabels );

options.verbose = 0;
if options.verbose
    for t=1:data.nframe
        map(:,:,t) = s_fpro(data.superpixels{t});        
    end
end
all_locains = sub2ind([data.height data.width], ceil(data.centres(:,1)), ceil(data.centres(:,2))); 
s_hog = zeros(data.slabels,9);
for t=1:data.nframe
    I=single(data.frames{t})/255;
    H=hog(I,8,9);
    H=imresize(H, [data.height data.width], 'nearest');
    H=reshape(H,[data.height*data.width,9]);
    loins = all_locains(data.bounds(t):data.bounds(t+1)-1);
    s_hog(data.bounds(t):data.bounds(t+1)-1,:) = H(loins,:);
end
s_lo = [(data.centres(:,1)-data.height/2)/(data.height/2) (data.centres(:,2)-data.width/2)/(data.width/2) ];
color_span = 256/options.bins;
[ s_colHist1, s_colHist2, s_colHist3 ] = getSuperpixelColHist( cellfun(@(x){floor(double(x)/color_span)}, data.frames), data.superpixels, data.slabels, options.bins);
s_colHist = [s_colHist1 s_colHist2 s_colHist3];
s_features = [s_hog s_lo s_colHist];
F =15;
M = 4*(2*F+1);
% W = sparse(data.bounds(data.nframe)-1,data.bounds(data.nframe)-1);
x_index = [];
y_index = [];
weight = [];
for t = 1:data.nframe
    query_index = [data.bounds(t):data.bounds(t+1)-1];
    test_index = [data.bounds(max(t-F,1)):data.bounds(min(data.nframe,t+F)+1)-1];
    Mx = 4*(min(data.nframe,t+F)-max(t-F,1));
    query_data = s_features(query_index, :);
    test_data = s_features(test_index, :);
    [sort_Ind,sort_Dis] = KDSearch(test_data',query_data',Mx,1);
    sort_Ind = sort_Ind';
    sort_Dis = sort_Dis';
    sort_Ind = test_index(sort_Ind);
    query_index = repmat(query_index',[Mx 1]);
    x_index = [x_index;query_index(:)];
    y_index = [y_index;sort_Ind(:)];
    weight = [weight;exp(-sort_Dis(:))];  
end
W=sparse(x_index,y_index,weight,data.bounds(data.nframe+1)-1,data.bounds(data.nframe+1)-1);
iD = sparse(1:data.bounds(data.nframe+1)-1,1:data.bounds(data.nframe+1)-1,1./sum(W,2));
P = iD*W;
v = double(s_fpro);


video_label(:) = 1;
for T = 1:10
   v = P*v; 
   video_l = v(data.Allsuperpixels);
   video_l(:,:,1) = double(data.gt);
   video_l(fore_tr_loind) = 1;
   video_l(back_tr_loind) = 0;
   v = double(getSuperpixelFprob( double(video_l), video_label ,data.superpixels, data.slabels ));
end

for t=1:data.nframe
    seg(:,:,t) = v(data.superpixels{t})>=0.5;
    imwrite(seg(:,:,t),[options.outputfolder '\' data.names{t} '.png']);
    [~,~,imgMarkup]=segoutput(im2double( data.frames{t}),double(seg(:,:,t)));
    imwrite(imgMarkup,[options.outputfolder '\' 'xx_' data.names{t} '.png']);
end
end