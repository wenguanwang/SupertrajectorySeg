function [tr] = linkFlowLDOF(data, options)
p = data.height;
q = data.width;

nf = data.nframe;
sample_step=options.sample_step;
margin=options.margin;

[ys,xs]=ndgrid(margin:sample_step:p-margin,margin:sample_step:q-margin);
ys=ys(:);
xs=xs(:);
sinds=sub2ind([p,q],ys,xs);%将坐标（ys，xs）在原【p，q】大小矩阵中，转化成下标大小。得到是：采样的坐标

X=zeros(nf,4*ceil(p*q/(sample_step^2))); %sample_ste^2=patch大小;/=一帧内patch个数
Y=zeros(nf,4*ceil(p*q/(sample_step^2)));

R=zeros(nf,4*ceil(p*q/(sample_step^2))); 
G=zeros(nf,4*ceil(p*q/(sample_step^2)));
B=zeros(nf,4*ceil(p*q/(sample_step^2))); 


Start=zeros(4*ceil(p*q/(sample_step^2)),1);
Prev_prob = zeros(4*ceil(p*q/(sample_step^2)),1);
Curr_prob = zeros(4*ceil(p*q/(sample_step^2)),1);

Ini_C=zeros(4*ceil(p*q/(sample_step^2)),3); 
im = reshape(data.frames{1},[data.height*data.width 3]);
Ini_C(1:length(ys),:) = im(sinds,:);

X(1,1:length(xs))=xs; %%所有点轨迹的起始坐标。
Y(1,1:length(ys))=ys;

Start(1:length(ys),1)=1;
Xtr=cell(nf,1);
Ytr=cell(nf,1);
Ttr=cell(nf,1);
Rtr=cell(nf,1);
Gtr=cell(nf,1);
Btr=cell(nf,1);
lens=cell(nf,1);

[Ys_ori,Xs_ori]=ndgrid([1:p],[1:q]);
xx = Xs_ori(:); %Represents cols
yy = Ys_ori(:); %Represents rows
pinds = sub2ind([p q], Ys_ori(:), Xs_ori(:));
% cur_idx = sub2ind([p q], orig_y(:), orig_x(:));  
% prev_prob = zeros(p,q);    
for t=1:data.nframe-1
    progress(sprintf('\t\t compute trajectory in frame'),t,data.nframe-1);
    I1 = double(data.frames{t});
    I2 = double(data.frames{t+1});
    
    to_imgR = double(I1(:,:,1));
    to_imgG = double(I1(:,:,2));
    to_imgB = double(I1(:,:,3));
      
    tracking_flow = data.fflow{t};
    reverse_flow  = data.bflow{t};
%     next_flow = readFlowFile([para.flow_dir 'Forward' get_image_name(imnames(t+1).name)  'LDOF.flo']);           
    tracking_flow_u = tracking_flow(:,:,1);
    tracking_flow_v = tracking_flow(:,:,2);
       
    Inds=find(Start>0);
    
    xs1=X(t,Inds);
    ys1=Y(t,Inds);
    cur_idx = sub2ind([p q], round(ys1), round(xs1))'; 
    
    R(t,Inds)=to_imgR(cur_idx); %%记录点轨迹颜色
    G(t,Inds)=to_imgG(cur_idx);
    B(t,Inds)=to_imgB(cur_idx);
    
    U=tracking_flow(:,:,1);
    V=tracking_flow(:,:,2);
    
    xs1_new=xs1'+U(cur_idx);%%得到下一帧的点轨迹坐标。
    ys1_new=ys1'+V(cur_idx);  
    xs1_new(xs1_new <= 1) = 1;
    ys1_new(ys1_new <= 1) = 1;
    xs1_new(xs1_new > q) = q;
    ys1_new(ys1_new > p) = p;      
    cur_idx_new = sub2ind([p q],round(ys1_new),round(xs1_new));
    
% %     U(inds_off)=-Inf;
% %     V(inds_off)=-Inf;
%     %bilinear interpolation
%     u=interp2(U,xs1,ys1);
%     v=interp2(V,xs1,ys1);
%     xs2=xs1+u;%%得到下一帧的点轨迹坐标。
%     ys2=ys1+v;  
%     xs2(xs2 <= 1) = 1;
%     ys2(ys2 <= 1) = 1;
%     xs2(xs2 > q) = q;
%     ys2(ys2 > p) = p;      
%     idx_new = sub2ind([p q],round(ys2),round(xs2));
     
%     u=interp2(U,xx,yy);
%     v=interp2(V,xx,yy);
%     xx_new=xx+u;%%得到下一帧的点轨迹坐标。
%     yy_new=yy+v;  
    xx_new=xx+U(:);%%得到下一帧的点轨迹坐标。
    yy_new=yy+V(:);  
    xx_new(xx_new <= 1) = 1;
    yy_new(yy_new <= 1) = 1;
    xx_new(xx_new > q) = q;
    yy_new(yy_new > p) = p;      
    pinds_new = sub2ind([p q],round(yy_new),round(xx_new));
  
    curr_cost = zeros(p,q);              
    %Penalizing for change in appearance 
    for img_dim=1:3
         to_img = double(I1(:,:,img_dim));
         from_img = double(I2(:,:,img_dim));
         curr_cost(pinds) = curr_cost(pinds) + (to_img(pinds)-from_img(pinds_new)).^2/256/256/3; 
         curr_cost(cur_idx) = curr_cost(cur_idx) + (Ini_C(Inds,img_dim)-from_img(cur_idx_new)).^2/256/256/3;          
    end       
        
    %Penalizing for occlusions
%     Xs_new=Xs_ori+tracking_flow(:,:,1);
%     Ys_new=Ys_ori+tracking_flow(:,:,2);    
%     flo_back_interp(:,:,1) = interp2(reverse_flow(:,:,1), Xs_new, Ys_new);
%     flo_back_interp(:,:,2) = interp2(reverse_flow(:,:,2), Xs_new, Ys_new);
%     discrepancy = tracking_flow+flo_back_interp;
%     discrepancy = discrepancy(:,:,1).^2 + discrepancy(:,:,2).^2;
% %     discrepancy = discrepancy ./ (tracking_flow(:,:,1).^2 + ...
% %         tracking_flow(:,:,2).^2+reverse_flow(:,:,1).^2 + reverse_flow(:,:,2).^2+1); 
%     discrepancy = discrepancy ./ (tracking_flow(:,:,1).^2 + ...
%         tracking_flow(:,:,2).^2+flo_back_interp(:,:,1).^2 + flo_back_interp(:,:,2).^2+1); 
%     curr_cost = curr_cost + discrepancy;  
    
    discrepancy = zeros(p,q); 
    t_flow_u = tracking_flow(:,:,1);
    t_flow_v = tracking_flow(:,:,2);
    r_flow_u = reverse_flow(:,:,1);  
    r_flow_v = reverse_flow(:,:,2);  
    discrepancy(pinds) = (t_flow_u(pinds) + r_flow_u(pinds_new)).^2+(t_flow_v(pinds) + r_flow_v(pinds_new)).^2;
    discrepancy(pinds) = discrepancy(pinds)./ (t_flow_u(pinds).^2 + ...
        t_flow_v(pinds).^2+r_flow_u(pinds_new).^2 + r_flow_v(pinds_new).^2+1);
    curr_cost = curr_cost + discrepancy;    
                
    Curr_prob(Inds) = Prev_prob(Inds)+ ((1-Prev_prob(Inds)).*(1-exp(-curr_cost(cur_idx)))); 
    cprob = Curr_prob(Inds);
    ON = cprob>0.5;

    Curr_prob(Curr_prob>0.5) = 0;
    Prev_prob = Curr_prob; 
 
    U(cur_idx(ON))=-Inf;
    V(cur_idx(ON))=-Inf;
    
    U(cur_idx(~ON))= tracking_flow_u(cur_idx(~ON));
    V(cur_idx(~ON))= tracking_flow_v(cur_idx(~ON));    
    
    %bilinear interpolation
    u=interp2(U,xs1,ys1);
    v=interp2(V,xs1,ys1);
    xs2=xs1+u;
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
    [Xtr{t},Ytr{t},Ttr{t},Rtr{t},Gtr{t},Btr{t},lens{t}]=terminate_at_t(inds_off,t,Start,X,Y,R,G,B);
    Start(inds_off)=0;
      
    %sample new points at t+1-> cover with 1 the area to be sampled and 0 the area to be left untouched
    map_occupied=ones(p,q);
    pixel_inds=sub2ind([p q],round(ys2),round(xs2));
    map_occupied(pixel_inds)=0;
    if sample_step>1
        map_occupied=imerode(map_occupied,strel('disk',sample_step));
        map_occupied=double(map_occupied>0);
    end
    ys_add=ys(map_occupied(sinds)>0);
    xs_add=xs(map_occupied(sinds)>0);
    
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
    im = reshape(data.frames{t+1},[data.height*data.width 3]);
    Ini_C(free_inds(1:num_add),:) = im(sub2ind([p,q],ys_add(:),xs_add(:)),:);
end
t=t+1;
%terminate
inds_on=find(Start>0);
[Xtr{t},Ytr{t},Ttr{t},Rtr{t},Gtr{t},Btr{t},lens{t}]=terminate_at_t(inds_on,t,Start,X,Y,R,G,B);
% [Xtr{t},Ytr{t},Ttr{t},lens{t}]=terminate_at_t(inds_on,t,Start,X,Y);

A=[Xtr{:};Ytr{:};Ttr{:};Rtr{:};Gtr{:};Btr{:}];
LL=cat(1,lens{:})';
B=mat2cell(A,6,LL);
tr=cell2struct(B,{'XYTRGB'},1);
tr=tr(LL>1);

% A=[];
% LL=cat(1,lens{:})';
% B=mat2cell(A,3,LL);
% trC=cell2struct(B,{'RGB'},1);
% trC=trC(LL>1);
end

function [Xt,Yt,Tt,Rt,Gt,Bt,lens_c]=terminate_at_t(inds_end,t,Start,X,Y,R,G,B)
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

Rt=R(inds);
Gt=G(inds);
Bt=B(inds);

Tt=ts;
end