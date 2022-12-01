function plot_connections_inactive_gray(f,img,rangein,csize,figno)
if ~isempty(figno)
    if nargin>2
        if ishandle(figno) close(figno); end
        figure(figno),
    end
end
x=f(1,:)';
y=f(2,:)';
nc=length(x);
weights=f(3,:)';
if ~isempty(img)
if size(img,3)>1
    imshow(im2double(rgb2gray(img)));
else
     imshow(img);
end
hold on;
end

S=csize*ones(nc,1);
cmap = colormap(jet);
colormap gray;
if isempty(rangein)
    if length(unique(weights))>1
        if min(weights)<0
            weights=weights+abs(min(weights));
        elseif min(weights)>0
            weights=weights-abs(min(weights));
        end
    end
    my_range=range(weights);
else
    my_range=rangein;

    if length(unique(weights))>1
        if min(weights)<0
            weights=weights+abs(min(weights));
        end
    end
end
if (my_range==0)
    weights=1;
    my_range=1;
end
cmap2 = cmap(floor((size(cmap,1)-1)*weights/my_range)+1,:);
scatter(x,y,S,cmap2,'filled');%weights','filled');
% hold on;
% axis image;
% hold on;
% axis image;
% colormap gray;
% hold off;

end
