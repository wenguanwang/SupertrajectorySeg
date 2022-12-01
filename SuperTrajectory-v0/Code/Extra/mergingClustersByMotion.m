
function geoNormDist = mergingClustersByMotion( Atr, tr, options )

    left = 1;
    right = 2;
    lLen = length( tr{ left } );
    tr = [ tr{ left }; tr{ right } ];
    
    [ trMat, trStart, trNum ] = formatTr2Mat( tr );
    
    % Compute the velocity of each trajectory    
%     [trVelocity patchMotion]= computeVelocity(tr, trLabels);
    trVelocity = computeVelocity( trMat, trStart, trNum );
    maxMotion = max( trVelocity( : ) );
    trMotion = trVelocity ./ maxMotion;

    % Find the boundary clusters
    avrCoordinates = computeAvrCoordinates( trMat, trStart, trNum );

    % find the top/bottom/left/right most label
    distThreshY = options.imgRow / 16;
    distThreshX = options.imgCol / 16;
    isMin = 1;
    topTrLabels = findBoundaryLabels( trMat, avrCoordinates, isMin, 2, distThreshY );
    leftTrLabels = findBoundaryLabels( trMat, avrCoordinates, isMin, 1, distThreshX );
    isMin = 0;
    bottomTrLabels = findBoundaryLabels( trMat, avrCoordinates, isMin, 2, distThreshY );
    rightTrLabels = findBoundaryLabels( trMat, avrCoordinates, isMin, 1, distThreshX );

    boundTrLabels = [ topTrLabels; bottomTrLabels; leftTrLabels; rightTrLabels ];
    boundTrLabels = unique( boundTrLabels );
    boundTrMotion = trMotion( boundTrLabels );
    
    
    % K-means Clustering
    % Set the clusters centres
    [ maxMotion, maxIdx ] = max( boundTrMotion );
    [ minMotion, minIdx ] = min( boundTrMotion );
    % K = 2
    centres = [maxIdx minIdx];
    [ centres, trIDX ] = kclusters( boundTrMotion', centres );

    % Detect the foreground patches in the detected background patches
    foreinbackLabels = detectForeinBackground( trIDX );

    % if there are some foreground patches, remove them from the background
    if ~isempty(foreinbackLabels)
        boundTrLabels(foreinbackLabels) = [];
    end

    % quoted on 2016/6/30
    frames = max(unique(trMat(3,:)));
    geoNormDist = computeGeoDist(tr, boundTrLabels, Atr, frames, para, lLen);

end

function trlabels = findBoundaryLabels( trMat, aveCoordinates, isMin, idx, thresh )
% find the boundary pathes' labels

    if( isMin )
        L = min( trMat( idx, : ) );
        Lstart = round( L + thresh ./ 2 );
        Lend = round( L + thresh );
    else
        L = max( trMat( idx, : ) );
        Lend = round( L - thresh ./ 2 );
        Lstart = round( L - thresh );
    end

    S = find( aveCoordinates( : , idx ) >= Lstart );
    E = find( aveCoordinates( : , idx ) <= Lend );
    labels = intersect( S, E );
    trlabels = unique( labels );

end

function trVelocity = computeVelocity( trMat, trStart, trNum )
% Compute the velocity of trajectories
    
    trVelocity = zeros( 1,trNum );
    for( i = 1 : trNum )
        trLen = trStart( 1, i + 1 ) - trStart( 1, i );
        corX = 0;
        corY = 0;

        % sum the total coordinates of the trajectory
        for( j = 1 : trLen-1 )
            index = trStart( 1, i ) + j;
            corX = corX + abs( trMat( 1, index ) - trMat( 1, index - 1 ) );
            corY = corY + abs( trMat( 2, index ) - trMat( 2, index - 1 ) );                
        end

        % displacement
        dspm = sqrt( corX ^ 2 + corY ^ 2 );
        trVelocity( 1, i ) = dspm ./ trLen;
        if( trVelocity( 1, i ) < 0.009 )
            trVelocity( 1, i ) = 0;
        end
        
    end

end

function avrCoordinates = computeAvrCoordinates( trMat, trStart, trNum )

    avrCoordinates = zeros( trNum, 2 );

    for( i = 1 : trNum )
        trLen = trStart( 1, i + 1 ) - trStart( 1, i );
        corX = 0;
        corY = 0;
        loopcnt = 0;
        % Sum the total coordinates of the trajectory
        for( framecnt = 0 : trLen - 1 )
            index = trStart( 1, i ) + framecnt;
            corX = corX + trMat( 1, index );
            corY = corY + trMat( 2, index );
            loopcnt = loopcnt + 1;
        end
        avrCoordinates( i, 1 ) = corX ./ loopcnt;
        avrCoordinates( i, 2 ) = corY ./ loopcnt;
    end

end

function foreLabels = detectForeinBackground( IDX )
% this function used to detect whether there are foreground patches in
% detected background patches
% input:
% IDX: the cluster result of the oversegmented patches

    clusterNum = unique( IDX );
    patchNum = length( IDX );

    if clusterNum ~= 2
        disp( 'The cluster Number must be 2!\n' );
        error(0);
    end

    foreLabels = [ ];
    label1Num = find( IDX == 1 );
    label2Num = find( IDX == 2 );
    foreinbackRatio1 = length( label1Num ) ./ patchNum;
    foreinbackRatio2 = length( label2Num ) ./ patchNum;
    if( foreinbackRatio1 <= 0.25 )
        foreLabels = label1Num;
    elseif( foreinbackRatio2 <= 0.25 )
        foreLabels = label2Num;
    end

end

function [ dictionary, index ] = kclusters( data, centres )
% This function clusters the input data into n sub-clusters
% n equals to the number of centres
% Input:
% data: to be clustered
% centres: clusters' centres

    centres = data( centres, : );

    featuretype = 'kclusters.mat';
    niters = 100;
    [ ncentres, dim ] = size( centres );
    old_centres = centres;
    display( 'Run k-means' );

    for( n = 1 : niters )
        % Save old centres to check for termination
        e2 = max( max( abs( centres - old_centres ) ) );

        inError( n ) = e2;
        old_centres = centres;
        tempc = zeros( ncentres, dim );
        num_points = zeros( 1, ncentres );

        id = eye( ncentres );
        d2 = EuclideanDistance( data, centres );
        % Assign each point to nearest centre
        [ minvals, index ] = min( d2', [ ], 1 );
        % matrix, if word i is in cluster j, post(i,j)=1, else 0;
        post = id( index, : ); 

        num_points = num_points + sum( post, 1 );

        for( j = 1 : ncentres )
            tempc( j, : ) =  tempc( j, : ) + sum( data( find( post( :, j ) ), : ), 1 );
        end

        for( j = 1 : ncentres )
            if( num_points( j ) > 0 )
                centres( j, : ) =  tempc( j, : ) / num_points( j );
            end
        end
        if n > 1
        % Test for termination
        
            %Threshold
            ThrError = 0.009;

            if( max( max( abs( centres - old_centres ) ) ) < 0.009 )
                dictionary = centres;
                fprintf('Saving texton dictionary\n');
                save ( [ './',featuretype ], 'dictionary' );      % save the settings of descriptor in opts.globaldatapath
                break;
            end
            fprintf( 'The %d th interation finished \n',n );
        end

    end
    
end

function geoNormDist = computeGeoDist( tr, boundLabels, Atr, frames, lLen )
%   Compute the geodesic saliency distance of each frame

    left = 1;
    right = 2;
    boundTr = tr( boundLabels );
    [ boundMat, boundStart, boundNum ] = formatTr2Mat( boundTr );
    [ trMat,trStart, trNum ] = formatTr2Mat( tr );
    
    boundMap = zeros( frames, boundNum );
    AtrMap = zeros( frames, trNum );
    
    geoNormDist = cell( 2, 1 );
    geoNormDist{ left } = cell( frames, 1 );
    geoNormDist{ right } = cell( frames, 1 );
       
    for i = 1 : boundNum
    %   Find the subset trajectories that crosses the current frame
        index = [ boundStart( 1, i ) : boundStart( 1, i + 1 ) - 1 ];
        frameid = boundMat( 3, index );
        boundMap( frameid, i ) = 1;
    end
    
    for i = 1 : trNum
    %   Find the subset trajectories that crosses the current frame 
        index = [ trStart( 1, i ) : trStart( 1, i + 1 ) - 1 ];
        frameid = trMat( 3, index );
        AtrMap( frameid, i ) = 1;
    end
    
    for i = 1:frames
    %   find the subset trajectories in each frame
        boundTrSub = boundLabels( boundMap( i, : ) > 0 );
        AtrSub = AtrMap( i, : )' * AtrMap( i, : );        
        
%         adjcMatrix = Atr>0;
        adjcMatrix = Atr .* AtrSub > 0;
        bdIds = boundTrSub;
%         colDistM = 1 - Atr;
        colDistM = 1 - Atr .* AtrSub;
        clip_value = 0;
        un = 0;
        
        geoDist = GeodesicSaliency( adjcMatrix, bdIds, colDistM, clip_value, un );
        % Delete the infinite value
        geoDist( find( geoDist == inf ) ) = 0;
        maxdist = max( geoDist( : ) );
        geoNormDist{ left }{ i } = geoDist( 1 : lLen ) ./ maxdist;
        geoNormDist{ right }{ i } = geoDist( lLen + 1 : end ) ./ maxdist;

    end    
    
end

function [foreNormDist] = computeForeGeoDist(tr, boundLabels, Atr, frames, para, lLen)
%   compute the geodesic saliency distance of each frame
    left = 1;
    right = 2;
%     boundTr = tr(boundLabels);
%     [boundMat boundStart boundNum] = formatTr2Mat(boundTr);
%     [trMat trStart trNum] = formatTr2Mat(tr);
    
%     boundMap = zeros(frames, boundNum);
%     AtrMap = zeros(frames, trNum);
    
    foreNormDist = cell(2, 1);
%     foreNormDist{left} = cell(frames, 1);
%     foreNormDist{right} = cell(frames, 1);
    
    
%     for i = 1:boundNum
%     %   find the subset trajectories that crosses the current frame
%         index = [boundStart(1,i): boundStart(1,i+1)-1];
%         frameid = boundMat(3,index);
%         boundMap(frameid,i) = 1;
%     end
%     
%     for i = 1:trNum
%     %   find the subset trajectories that crosses the current frame 
%         index = [trStart(1,i):trStart(1,i+1)-1];
%         frameid = trMat(3,index);
%         AtrMap(frameid,i) = 1;
%     end
    
%     for i = 1:frames
    %   find the subset trajectories in each frame
%         boundTrSub = boundLabels(boundMap(i,:)>0);
%         AtrSub = AtrMap(i,:)'*AtrMap(i,:);
        
        
        adjcMatrix = Atr>0;
%         adjcMatrix = Atr.*AtrSub>0;
        bdIds = boundLabels;
        colDistM = 1 - Atr;
%         colDistM = 1 - Atr.*AtrSub;
        clip_value = 0;
        un = 0;
        
        geoDist = GeodesicSaliency(adjcMatrix, bdIds, colDistM, clip_value, un);
        % delete the infinite value
        geoDist(find(geoDist == inf)) = 0;
        maxdist = max(geoDist(:));
%         geoNormDist{i} = geoDist./maxdist;
        foreNormDist{left} = 1 - geoDist(1:lLen)./maxdist;
        foreNormDist{right} = 1 - geoDist(lLen+1:end)./maxdist;
%         geoDist = [];
%     end
    
    
end

