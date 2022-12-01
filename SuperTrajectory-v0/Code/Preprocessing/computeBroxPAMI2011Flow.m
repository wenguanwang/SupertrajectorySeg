% Wrapper to compute the broxPAMI2011 optical flow of a given shot
%
%    Copyright (C) 2013  Anestis Papazoglou
%
%    You can redistribute and/or modify this software for non-commercial use
%    under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    For commercial use, contact the author for licensing options.
%
%    Contact: a.papazoglou@sms.ed.ac.uk

function [ flow, backflow ] = computeBroxPAMI2011Flow( options, shot )
    
    startFrame = options.ranges( shot );
    endFrame = options.ranges( shot + 1 ) - 1;
    

    
    [ flow{ left }, backflow{ left } ] = flowGeneration( options, startFrame, endFrame, left );


end

function [ flow, backflow ] = flowGeneration( options, startFrame, endFrame, stereo )
%   generate the flow of stereo videos

    flowframes = endFrame - startFrame;
    flow = cell( 1, flowframes );
    backflow = cell( 1, flowframes );
    
    totalTimeTaken = 0;
    
    currImage = readFrame( options, startFrame );
    if( size( currImage, 3 ) == 1 )
        currImage = gray2rgb( currImage );
    end
    currImage = double( currImage );
    for( i = startFrame + 1: endFrame )
        tic
        
        index = i - startFrame;
        nextImage = readFrame( options, i );
        if( size( nextImage, 3 ) == 1 )
            nextImage = gray2rgb( nextImage );
        end
        nextImage = double( nextImage );

        if( options.vocal )
            fprintf( 'computeBroxPAMI2011Flow: Computing optical flow of pair: %i of %i... ', ...
                index, flowframes );
        end

        flowframe = mex_LDOF( currImage, nextImage );
        backflowframe = mex_LDOF( nextImage, currImage );
        
        flowmap = flowToColor(flowframe);
        
%         flowframe = int16( flowframe );
%         backflowframe = int16( backflowframe );
        
        flowframe = single( flowframe );
        backflowframe = single( backflowframe );
        
        
        flowmap = flowToColor(double(flowframe));
        
        flow{ index }( :, :, 1 ) = flowframe( :, :, 2 );
        flow{ index }( :, :, 2 ) = flowframe( :, :, 1 );
        backflow{ index }( :, :, 1 ) = backflowframe( :, :, 2 );
        backflow{ index }( :, :, 2 ) = backflowframe( :, :, 1 );
        
        
        currImage = nextImage;
        timeTaken = toc;
        totalTimeTaken = totalTimeTaken + timeTaken;
        
        if( options.vocal )
            fprintf( 'done. Time taken: %.2f sec\n', timeTaken );
        end
        
    end
    
    if( options.vocal )
        if( options.isLeft )
            reminder = 'computeBroxPAMI2011Flow: (Left Flow) Total time taken: %.2f sec\n';
        else
            reminder = 'computeBroxPAMI2011Flow: (Right Flow) Total time taken: %.2f sec\n';
        end
        fprintf( reminder, totalTimeTaken );
        fprintf( 'computeBroxPAMI2011Flow: Average time taken per frame: %.2f sec\n', ...
            totalTimeTaken / flowframes );
    end

end

function flow = flowBackwardGeneration( options, startFrame, endFrame, stereo )
%   generate the flow of stereo videos
    flowframes = endFrame - startFrame;
    flow = cell( 1, flowframes );
    
    totalTimeTaken = 0;
    
    currImage = readFrame( options, startFrame );
    if( size( currImage, 3 ) == 1 )
        currImage = gray2rgb( currImage );
    end
    currImage = double( currImage );
    for( i = startFrame + 1: endFrame )
        tic
        
        index = i - startFrame;
        nextImage = readFrame( options, i );
        if( size( nextImage, 3 ) == 1 )
            nextImage = gray2rgb( nextImage );
        end
        nextImage = double( nextImage );

        if( options.vocal )
            fprintf( 'computeBroxPAMI2011Flow: Computing optical flow of pair: %i of %i... ', ...
                index, flowframes );
        end

        flowbackframe = mex_LDOF( nextImage, currImage );
        flowmap = flowToColor(flowbackframe);
        flowbackframe = int16( flowbackframe );
        flowmap = flowToColor(double(flowbackframe));
        flow{ index }( :, :, 1 ) = flowbackframe( :, :, 2 );
        flow{ index }( :, :, 2 ) = flowbackframe( :, :, 1 );

        currImage = nextImage;
        timeTaken = toc;
        totalTimeTaken = totalTimeTaken + timeTaken;
        
        if( options.vocal )
            fprintf( 'done. Time taken: %.2f sec\n', timeTaken );
        end
        
    end
    
    if( options.vocal )
        if( options.isLeft )
            reminder = 'computeBroxPAMI2011Flow: (Left Backward Flow) Total time taken: %.2f sec\n';
        else
            reminder = 'computeBroxPAMI2011Flow: (Right Backward Flow) Total time taken: %.2f sec\n';
        end
        fprintf( reminder, totalTimeTaken );
        fprintf( 'computeBroxPAMI2011Flow: (Backward) Average time taken per frame: %.2f sec\n', ...
            totalTimeTaken / flowframes );
    end

end
