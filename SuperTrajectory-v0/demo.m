clc
clear all
warning off;
addpath( genpath( '.' ) );
options.sample_step = 1;
options.margin = 1;
options.aggr=5;                              %% velocity support
options.my_var_euk=0.01;                                       %% euclidean variance
options.my_var_cuk=0.1;                                       %% euclidean variance
options.my_var_ut=0.1; 
options.str_num=1000;

options.regnum = 2000;
options.m = 10;
options.bins = 20;

options.verbose = 1;
options.m1 = 1;
options.m2 = 1;
options.foldername = fileparts( mfilename( 'fullpath' ) );
options.datasetname = 'inputs';
videoFiles = dir(fullfile(options.foldername, 'Data', options.datasetname));
videoNUM = length(videoFiles)-2;


for videonum = 1:videoNUM
    
    folder =  videoFiles(videonum+2).name
    options.infolder = fullfile(options.foldername, 'Data', options.datasetname, folder);
    options.datafolder = fullfile(options.foldername, 'Data', options.datasetname, folder,'data');
    options.outputfolder = fullfile(options.foldername, 'Data', options.datasetname, folder,'result');
       
    if( ~exist( options.datafolder, 'dir' ) ) 
         mkdir( options.datafolder ); 
    end;
    if( ~exist( options.outputfolder, 'dir' ) ) 
         mkdir( options.outputfolder ); 
    end;
    
    % Cache all frames in memory
    [data.frames,data.names,data.height,data.width,data.nframe ]= readAllFrames( options.infolder );   
    data.gt = imread(fullfile(options.foldername, 'Data','Annotations', folder,[data.names{1} '.png']));
    data.gt = imresize(data.gt, [data.height data.width])>122;
    
    %     [data.gtframes,data.gtnames,~,~,~ ]= readAllFrames( fullfile(options.foldername, options.datasetname, videoFiles(videonum+2).name, 'A') );
    [superpixelsLAB, data.superpixelsLabel]= loadSuperpixels( options );
    if(isempty(data.superpixelsLabel ) )
        [superpixelsLAB, data.superpixelsLabel]= computeSuperpixels( options, data.frames );
    end  
    [ data.superpixels, data.bounds, data.slabels ] = makeSuperpixelIndexUnique( data.superpixelsLabel );     
    [ colours, data.centres, ~ ] = getSuperpixelStats( cellfun(@(x){double(x)}, data.frames), data.superpixels, data.slabels );
    
    x = data.centres(:,1);
    x(x<=0)=1;
    x(x>data.height)=data.height;
    y = data.centres(:,2);
    y(y<=0)=1;
    y(y>data.width)=data.width;
    
    data.centres = [x y];
    
    [data.fflow, data.bflow] = computeFlowLDOF(data, options);
    data.Allsuperpixels = zeros(data.height,data.width,data.nframe);
    data.AllframesR = zeros(data.height,data.width,data.nframe);
    data.AllframesG = zeros(data.height,data.width,data.nframe);
    data.AllframesB = zeros(data.height,data.width,data.nframe);
    for i = 1:data.nframe
        data.Allsuperpixels(:,:,i) = data.superpixels{i};
        data.AllframesR(:,:,i) = data.frames{i}(:,:,1);
        data.AllframesG(:,:,i) = data.frames{i}(:,:,2);
        data.AllframesB(:,:,i) = data.frames{i}(:,:,3);
    end
    data.Allfflowx = zeros(data.height,data.width,data.nframe-1);
    data.Allfflowy = zeros(data.height,data.width,data.nframe-1);
    data.Allbflowx = zeros(data.height,data.width,data.nframe-1);
    data.Allbflowy = zeros(data.height,data.width,data.nframe-1);
    for i = 1:data.nframe-1           
        data.Allfflowx(:,:,i) = data.fflow{i}(:,:,1);
        data.Allfflowy(:,:,i) = data.fflow{i}(:,:,2);
        data.Allbflowx(:,:,i) = data.bflow{i}(:,:,1);
        data.Allbflowy(:,:,i) = data.bflow{i}(:,:,2);
    end  
    
%% super point trajectory
    [tr, tr_lo, tr_mo, tr_co, all_center,all_cl]=super_point_trajectory(data, options);
    
 %% Compute Super-Trajectory Feature
    disp('compute super-trajectory feature');
    [str_tr, str_co, str_mo, str_lo] =get_str_Feature(tr_lo, tr_mo, tr_co, all_cl);
    options.verbose = 0;
    if options.verbose   
        t = 1;      
        [tr_XYT,tr_id]=quick_tr(tr);
        tr_ins = find(tr_XYT(3,:)==t);
        str_ins = all_cl(tr_id(tr_ins));
        tr_lo = tr_XYT(1:2,tr_ins)';
        or_img = zeros([data.height*data.width,3]);
        tr_loind = sub2ind([data.height data.width], int16(tr_lo(:,2)), int16(tr_lo(:,1)));
        or_img(tr_loind,:)=str_co(str_ins,:);
        or_img = uint8(reshape(or_img,[data.height,data.width,3]));
        imshow(or_img);
    end   
    
%% super point trajectory segmentation
    [seg]=seg_str(data, options,tr,all_cl,str_tr, str_lo, str_mo, str_co);
end
