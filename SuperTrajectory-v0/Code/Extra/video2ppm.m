function [video_dir left_video_dir right_video_dir] = video2ppm(orivideo_dir, video_name)
% input:
% orivideo_dir: is the directory of original video
% output:
% video_dir: is the processed video directory

% read the original video
video = VideoReader(orivideo_dir);
video_n = video_name(1:end-4);
frame_number = floor(video.Duration * video.FrameRate);

% video_dir = [pwd '\extra\video2jpg\' video_n '\'];
video_dir = [pwd '\extra\video2ppm\' video_n '\'];
if ~exist(video_dir,'dir') 
    mkdir(video_dir);
end 

% get the absolute size of the input video
[top bottom col] = get_video_size(orivideo_dir);

%% Split Frames into PPM Format
image_dir_l = [video_dir 'left\'];
image_dir_r = [video_dir 'right\'];
if ~exist(image_dir_l,'dir')
    mkdir(image_dir_l);
end
if ~exist(image_dir_r,'dir')
    mkdir(image_dir_r);
end
for i = 1 : 8 %8 :frame_number
    image_name = [video_dir video_n num2str(i)];
    image_name = strcat(image_name,'.ppm');
    
    image_name_l = [image_dir_l video_n '_left_' num2str(i)]; % left name
    image_name_r = [image_dir_r video_n '_right_' num2str(i)]; % right name
    image_name_l = strcat(image_name_l, '.ppm');
    image_name_r = strcat(image_name_r, '.ppm');
    
    I = read(video,i);                               % read the original frame
    I = I(top:bottom, :, :);                          % cut the frame into absolute
    I_left = I(:, 1:col,:);
    I_right = I(:, col+1:col*2,:);
    imwrite(I,image_name,'ppm');                   % write the absolute left and right frame
    imwrite(I_left,image_name_l,'ppm');                   
    imwrite(I_right,image_name_r,'ppm');                   
    clear I I_left I_right;
end
left_video_dir = [video_dir 'left\'];
right_video_dir = [video_dir 'right\'];
end

