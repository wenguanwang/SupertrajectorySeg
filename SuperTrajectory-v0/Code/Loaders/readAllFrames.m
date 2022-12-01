function  [frames, names, height,width,nframe ] = readAllFrames(infolder)

    frameFiles = imdir(fullfile(infolder));
    nframe = length( frameFiles);
    frames = cell( nframe, 1 );
    names = cell( nframe, 1 );
    for index = 1: nframe 
        [~, frameName] = fileparts(frameFiles(index).name);
        if exist(fullfile(infolder, [frameName '.png']),'file')
            frame = imread(fullfile(infolder, [frameName '.png']));
        elseif exist(fullfile(infolder, [frameName '.jpg']),'file')
            frame = imread(fullfile(infolder, [frameName '.jpg']));
        elseif exist(fullfile(infolder, [frameName '.bmp']),'file')
            frame = imread(fullfile(infolder, [frameName '.bmp']));
        end
        frames{ index } = frame;
        names{ index } =  frameName;
    end
    [ height,width ] = size(frame(:,:,1));
end
