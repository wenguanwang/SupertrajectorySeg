function tr=linkFlow(data, options)

p = data.height;
q = data.width;
[orig_y,orig_x] = meshgrid(1:p, 1:q);
nf = data.nframe;

% sample_step=options.sample_step;
% margin=options.margin;

% [ys,xs]=ndgrid(margin:sample_step:p-margin,margin:sample_step:q-margin);
[ys,xs]=ndgrid(1:p,1:q);
ys=ys(:);
xs=xs(:);

pinds=sub2ind([p,q],ys,xs);%将坐标（ys，xs）在原【p，q】大小矩阵中，转化成下标大小。得到是：采样的坐标
% X=zeros(nf,4*ceil(p*q/(sample_step^2))); %sample_ste^2=patch大小;/=一帧内patch个数
% Y=zeros(nf,4*ceil(p*q/(sample_step^2)));
% Start=zeros(4*ceil(p*q/(sample_step^2)),1);

X(1,1:length(xs))=xs; %%所有点轨迹的起始坐标。
Y(1,1:length(ys))=ys;
Start(1:length(ys),1)=1;

Xtr=cell(nf);
Ytr=cell(nf);
Ttr=cell(nf);
lens=cell(nf);

xx = orig_x(:); %Represents cols
yy = orig_y(:); %Represents rows

idx = sub2ind([p q], orig_y(:), orig_x(:));   
[Ys_ori,Xs_ori]=ndgrid([1:p],[1:q]);
prev_prob = zeros(p,q);    


for t=1:data.nframe-2
    progress(sprintf('\t\t compute trajectory in frame'),t,data.nframe-1);
    I1 = double(data.frames{t});
    I2 = double(data.frames{t+1});
    Ff = data.fflow{t};
    Fb = data.bflow{t};
    tracking_flow = Ff;
    reverse_flow  = Fb;
%     next_flow = readFlowFile([para.flow_dir 'Forward' get_image_name(imnames(t+1).name)  'LDOF.flo']);           
 
    tracking_flow_u = tracking_flow(:,:,1);
    tracking_flow_v = tracking_flow(:,:,2);
        
    xx_new = round(xx + tracking_flow_u(idx));
    yy_new = round(yy + tracking_flow_v(idx));
        
    curr_cost = zeros(p,q);
                
    xx_new(xx_new <= 0) = 1;
    yy_new(yy_new <= 0) = 1;
    xx_new(xx_new > q) = q;
    yy_new(yy_new > p) = p;  
        
    idx_new = sub2ind([p q], yy_new, xx_new);    
                          
    %Penalizing for change in appearance        
    for img_dim=1:3
         to_img = double(I1(:,:,img_dim));
         from_img = double(I2(:,:,img_dim));
         curr_cost(idx) = curr_cost(idx) + (to_img(idx)-from_img(idx_new)).^2/256/256/3;                
    end       
        
%     for flow_dim=1:2
%         to_flow = double(tracking_flow(:,:,flow_dim));
%         from_flow = double(next_flow(:,:,flow_dim));            
%         val = 1+to_flow(idx).^2+from_flow(idx_new).^2;                
%         curr_cost(idx) = curr_cost(idx) + (to_flow(idx)-from_flow(idx_new)).^2./val*5;        
%     end

    %Penalizing for occlusions
    Xs_new=Xs_ori+tracking_flow(:,:,1);
    Ys_new=Ys_ori+tracking_flow(:,:,2);    
    flo_back_interp(:,:,1) = interp2(reverse_flow(:,:,1), Xs_new, Ys_new);
    flo_back_interp(:,:,2) = interp2(reverse_flow(:,:,2), Xs_new, Ys_new);
    

%     flo_back_interp(:,:,1) = interp2(reverse_flow(:,:,1), xx_new, yy_new);
%     flo_back_interp(:,:,2) = interp2(reverse_flow(:,:,2), xx_new, yy_new);

% 	discrepancy = tracking_flow + reverse_flow;
    discrepancy = tracking_flow+flo_back_interp;
    discrepancy = discrepancy(:,:,1).^2 + discrepancy(:,:,2).^2;
    discrepancy = discrepancy ./ (tracking_flow(:,:,1).^2 + ...
        tracking_flow(:,:,2).^2+reverse_flow(:,:,1).^2 + reverse_flow(:,:,2).^2+0.0001); 
    curr_cost = curr_cost + discrepancy;    
                
    curr_prob = prev_prob + ((1-prev_prob).*(1-exp(-curr_cost)));                
        
        
    xx = xx_new;
    yy = yy_new;
    prev_prob = curr_prob;
%     to_frame = from_frame;
    
    ON = curr_prob<=0.5;
    
%     [ON]=prune_flow_field(I1,I2,Ff,Fb,para,0);%%去除光流的区域。
    inds_off=find(ON==0);
%     length(inds_off);
    U=Ff(:,:,1);
    V=Ff(:,:,2);
    U(inds_off)=-Inf;
    V(inds_off)=-Inf;
    prev_prob(inds_off) = 0;
    %bilinear interpolation
    Inds=find(Start>0);
    xs1=X(t,Inds);
    ys1=Y(t,Inds);
    u=interp2(U,xs1,ys1);
    v=interp2(V,xs1,ys1);
    xs2=xs1+u;%%得到下一帧的点轨迹坐标。
    ys2=ys1+v;
    
    
    %filter pts based on flow filtering
    choose=xs2>margin & ys2>margin & xs2<q-margin & ys2<p-margin;%%去除一出边界的点轨迹。
    inds_on=Inds(choose);
    xs2=xs2(choose);
    ys2=ys2(choose);
    X(t+1,inds_on)=xs2; %%将下一帧的位置加进去。
    Y(t+1,inds_on)=ys2;
    
    %terminate
    inds_off=setdiff(Inds,inds_on);
    [Xtr{t},Ytr{t},Ttr{t},lens{t}]=terminate_at_t(inds_off,t,Start,X,Y);
    Start(inds_off)=0;
    
    
    %sample new points at t+1-> cover with 1 the area to be sampled and 0 the area to be left untouched
    map_occupied=ones(p,q);
    pixel_inds=sub2ind([p q],round(ys2),round(xs2));
    map_occupied(pixel_inds)=0;
    if sample_step>1
        map_occupied=imerode(map_occupied,strel('disk',sample_step));
        map_occupied=double(map_occupied>0);
    end
    ys_add=ys(map_occupied(pinds)>0);
    xs_add=xs(map_occupied(pinds)>0);
    
%     if 0
%         figure(12),
%         imagesc(map_occupied)
%         hold on;
%         plot(xs_add,ys_add,'.g','markersize',0.1);
% %         plot(xs2,ys2,'.r');
%         axis image
%         title('newly added points')
%     end
    
    num_add=length(xs_add);
    free_inds=find(Start==0);
    X(t+1,free_inds(1:num_add))=xs_add;
    Y(t+1,free_inds(1:num_add))=ys_add;
    Start(free_inds(1:num_add))=t+1;
end
t=t+1;
%terminate
inds_on=find(Start>0);
[Xtr{t},Ytr{t},Ttr{t},lens{t}]=terminate_at_t(inds_on,t,Start,X,Y);

A=[Xtr{:};Ytr{:};Ttr{:}];
LL=cat(1,lens{:})';
B=mat2cell(A,3,LL);
tr=cell2struct(B,{'XYTPos'},1);

tr=tr(LL>1);

end

function [Xt,Yt,Tt,lens_c]=terminate_at_t(inds_end,t,Start,X,Y)
inds1=sub2ind(size(X),Start(inds_end),inds_end);
inds2=sub2ind(size(X),t*ones(length(inds_end),1),inds_end);
lens_c=t-Start(inds_end)+1;
L=[0;cumsum(lens_c)];
inds=zeros(1,sum(lens_c));
ts=zeros(1,sum(lens_c));
for i=1:length(inds1)
    inds(L(i)+1:L(i+1))=inds1(i):inds2(i);
    ts(L(i)+1:L(i+1))=Start(inds_end(i)):t;
end
Xt=X(inds);
Yt=Y(inds);
Tt=ts;
end

