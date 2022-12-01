function visualize_trajectory_affinity_graph(tr,V,...
    data,A,R,range_init,normalized,figno,id)

if isempty(V)
    V=ones(length(tr),1);
end
if ~exist('id','var')
    id=[];
end
[XYT,tr_id]=quick_tr(tr);
Lengths=hist(tr_id,[1:length(tr)]);

nf=data.nframe;
if ishandle(figno) close(figno); end
figure(figno),

ppos =spp(ceil(sqrt(nf)),ceil(sqrt(nf)));
minimum=min(V);
for ii=1:data.nframe
    img=data.frames{ii};
    h=subplot('position',ppos(ii,:));
    t=ii;
    %   h=subplot2(ceil(sqrt(nf)),ceil(sqrt(nf)),counter);
    inds=find(XYT(3,:)==t);
    x=XYT(1,inds);
    y=XYT(2,inds);
    tr_id_c=tr_id(inds);
    
    if minimum<0
        weightsc=V(tr_id_c)'+abs(minimum);
    else
        weightsc=V(tr_id_c)'-minimum;
    end
    plot_connections_inactive_gray([x; y; weightsc],img,range_init,15,[]);
    title(['frame:' num2str(t)]);
    hold on;
    
end

info.tr_ids_clicked=[];
info.XYT=XYT;
info.tr_id=tr_id;
info.pts_clicked=[];
info.ppos=ppos;
info.numhs=[];
set(gcf,'UserData',info);
set(gcf, 'WindowButtonDownFcn', @callback_click);
hold on;

cmap = colormap(jet);
colormap(gray)
fprintf('click on trajectory points to see their motion affinities..')

    function callback_click(src, event);
        info=get(src,'UserData');
        myaxis = get(gca,'position');
        [dum, subplotid] = min(sum(abs(ppos - repmat(myaxis,size(ppos,1),1)),2));
        pt = get(gca, 'CurrentPoint');
        %show clicked point
        pt=pt(1,:);
        %plot(pt(1),pt(2),'*y');
        hold on;
        for jj=1:length(info.numhs)
            delete(info.numhs(jj));
        end
        info.numhs=[];
        
        img=data.frames{subplotid};
%         imread(imnames(subplotid).name);
        t=subplotid;
        inds=find(XYT(3,:)==t);
        x=XYT(1,inds);
        y=XYT(2,inds);
        tr_id_c=tr_id(inds);
        
        %plot(x,y,'*g');
        
        f=[x;y;tr_id_c];
        pts=f';
        [minimum,minindex]=min((pts(:,1)-pt(1,1)).^2+(pts(:,2)-pt(1,2)).^2);
        pt=pts(minindex,1:2);
        plot(pt(1),pt(2),'*k');
        hold on;
        %show trajectory clicked
        string=sprintf('%d',pts(minindex,3));
        ht=text(pt(1),pt(2),string,'BackgroundColor',[.7 .9 .7]);
        info.numhs=[info.numhs;ht];
        if ~isempty(id)
            tr_id_cl=id;
        else
            tr_id_cl=f(3,minindex);
        end
        
        if ishandle(10) close(10); end
        figure(10),
        h2=subplot(3,1,1)
        plot_connections_inactive_gray([x;y;V(tr_id_c)'+1],img,[],10,[]);
        hold on;
        plot(pt(1),pt(2),'py','markersize',15);
        title(['frame ' num2str(t) 'length traj=' num2str(Lengths(tr_id_cl))]);
        if ~isempty(A)
            a=A(tr_id_cl,:);
            tr_ids_conn=find(a);
            [tr_ids_conn,ik1,ik2]=intersect(tr_id_c,tr_ids_conn);
            weights_attr=a(tr_ids_conn);
            fattr=[x(ik1);y(ik1);weights_attr];
            h3=subplot2(3,1,2);
            if size(tr_ids_conn,2)<size(tr_ids_conn,1)
                tr_ids_conn=tr_ids_conn';
            end
            plot_points_active_general([fattr; tr_ids_conn] ,img,222);
            figure(10),
            if normalized
                plot_connections_inactive_gray(fattr,img,1,20,[]);
            else
                plot_connections_inactive_gray(fattr,img,[],20,[]);
            end
            hold on;
            %axis image;
            plot(pt(1),pt(2),'py','markersize',15);
            %title('attractions');
        end
        
        if ~isempty(R)
            r=R(tr_id_cl,:);
            tr_ids_conn=find(r);
            [tr_ids_conn,ik1,ik2]=intersect(tr_id_c,tr_ids_conn);
            weights_rep=r(tr_ids_conn);
            frep=[x(ik1);y(ik1);weights_rep];
            subplot2(3,1,3);
            if normalized
                plot_connections_inactive_gray(frep,img,[],20,[]);
            else
                plot_connections_inactive_gray(frep,img,1,20,[]);
            end
            hold on;
            plot(pt(1),pt(2),'py','markersize',15);
           % title('repulsions');
        end
        hold on;
        subplot('position', [0.0,0.7,0.05,0.25]);
        colorbarimage = zeros(256,50,3);
        colorbarimage(:,:,1) = repmat(flipud(cmap(:,1)),1,50);
        colorbarimage(:,:,2) = repmat(flipud(cmap(:,2)),1,50);
        colorbarimage(:,:,3) = repmat(flipud(cmap(:,3)),1,50);
        imagesc(colorbarimage); axis image; axis off;
        %  strings={'-1' '1'};
        %  imagesc(colorbarimage); axis image; axis off;
        %  hold on;
        %  text(100,50,strings{2});
        %  hold on
        %  text(100,250,strings{1});
        set(src,'UserData',info);
    end

end