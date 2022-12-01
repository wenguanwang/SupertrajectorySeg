% Function to compute some given superpixel method for a given shot
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

function [superpixelsLAB,superpixelsLabel]= computeSuperpixels( options, frames )
    
    fprintf( 'computeSuperpixels: \n');
    filename = fullfile( options.datafolder, 'superpixels.mat' );     
    nframes =length(frames);
        
    totalTimeTaken = 0;

    superpixelsLAB = cell( nframes, 1 );
    superpixelsLabel = cell( nframes, 1 );
    [ height,width ] = size(frames{1}(:,:,1));
    PixNum = height*width;    
    for index = 1: nframes
            tic;
            fprintf( 'computeSLIC: Processing frame %i/%i... ', ...
                index, nframes );          
            frame = frames{index};            
            frameVecR = reshape( frame(:,:,1)', PixNum, 1);
            frameVecG = reshape( frame(:,:,2)', PixNum, 1);
            frameVecB = reshape( frame(:,:,3)', PixNum, 1); 
            frameAttr=[ height ,width, options.regnum, options.m, PixNum ];
            [ Label, Sup1, Sup2, Sup3, k ] = SLIC( double(frameVecR), double(frameVecG), double(frameVecB), frameAttr );
            Label = int32(reshape(Label,width,height)');
%             superpixels{ index }.Label = Label+1;
            superpixelsLAB{ index } = [Sup1 Sup2 Sup3];
            superpixelsLabel{ index } = Label+1;
            timeTaken = toc;
            totalTimeTaken = totalTimeTaken + timeTaken;
            fprintf( 'time taken: %.2f seconds\n', timeTaken );               
    end
    save( filename, 'superpixelsLAB', 'superpixelsLabel','-v7.3' );
    fprintf( 'computeSLIC: Total time taken: %.2f sec\n', totalTimeTaken );
    fprintf( 'computeSLIC: Average time taken per frame: %.2f sec\n', ...
                totalTimeTaken / nframes );  
    fprintf( 'computeSuperpixels: finished processing\n' );
end
