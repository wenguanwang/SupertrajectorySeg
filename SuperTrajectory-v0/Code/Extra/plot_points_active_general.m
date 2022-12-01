function plot_points_active_general(f,img,figno)
%herer we show the corresponding index from the 3rd row of f
%finding the indices of feature locations
%f is 2 X N:feature locations
line_inds=[];
%plots the lines in the cell lines and shows the number when you cklick on them
%and returns the numbers
%if strcmp(col,'rand')
if nargin>2
    if ishandle(figno) close(figno); end
    figure(figno),
end
imshow(img)
hold on;
plot(f(1,:),f(2,:),'.g','markersize',10);
%  for jj=1:size(f,2)
%  plot(f(1,jj),f(2,jj),'*g')
%  hold on;
%  end
hold on;
% Check whether we need to expand bounds


h=gcf;
st.f=f;
st.hs=[];
set(gca,'UserData',st);
set(h, 'WindowButtonDownFcn', @callback_click);
%line_inds_clicked=get(gca,'UserData');
    function callback_click(src, event);
        
        pt = get(gca, 'CurrentPoint');
        pt=pt(1,:);
        plot(pt(1),pt(2),'*b');
        hold on;
        selectionType = get(gcf, 'SelectionType');
        st=get(gca,'UserData');
        f=st.f;
        
        hs=st.hs;
        for ii=1:length(hs)
            delete(hs(ii));
        end
        
        hs=[];
        pts=f';
        [minimum,minindex]=min((pts(:,1)-pt(1,1)).^2+(pts(:,2)-pt(1,2)).^2);
        pt=pts(minindex,1:2);
        plot(pt(1),pt(2),'*k');
        hold on;
        pts=full(pts);
        string=sprintf('%0.3f',full((pts(minindex,3))));
        if size(pts,2)==4
            string=sprintf('%s \n %d',string,pts(minindex,4));
        end
        pt=double(pt);
        h=text(pt(1),pt(2),string,'BackgroundColor',[.7 .9 .7]);
        hs=[hs h];
        % For model pt selection
        %  if (strcmp(selectionType, 'normal'))
        %       for i=1:size(edgelist,2)
        
        %              if ismember(pt,edgelist{i},'rows')
        %                  i
        %                  %i found corresponding line segment
        %                  a=edgelist{i};
        %                  hold on;
        %                  plot(a(:,1),a(:,2),'gs','MarkerSize', 3);
        %                  hold on;
        %                  text(a(1,1),a(1,2),sprintf('%d',i),'BackgroundColor',[.7 .9 .7]);
        %                  line_inds_clicked=get(gca,'UserData');
        %                  line_inds=[line_inds_clicked;i];
        %  if ishandle(figno2) close(figno2); end
        %  figure(figno2),
        %  plot(features(:,i));
        %  hold off;
        %                  set(gca,'UserData',line_inds);
        %                  break;
        %              end
        %  end
        % end
        
        st.hs=hs;
        set(gca,'UserData',st);
    end

%   h2 = plot(x(pixel_order{contour_id}), y(pixel_order{contour_id}), 'gs', 'MarkerSize', 6);
%             h1 = plot(sx, sy, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
%             %             title(sprintf('Frame %d: eig %d, loop id=%d, [x,y]=[%d,%d]', ...
%             %                 frame_id, jj, loop_id, sx, sy));
%
%             title(sprintf('Frame %d: contour id=%d, [x,y]=[%d,%d]', ...
%                 frame_id, contour_id, sx, sy));

end
