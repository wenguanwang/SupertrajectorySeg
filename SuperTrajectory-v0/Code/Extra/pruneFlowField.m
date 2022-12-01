% Wrapper to compute the SLIC superpixels of a given shot
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

function [ ON ] = pruneFlowField( I1, Ff, Fb, debug )

    % Texture threshold1
    structurePortion = 0;
    % Texture threshold2
    minValue = 0.000005;
    % ForBackCheck 1
    fb_check = 0.1;
    % ForBackCheck 2
    fb_check2 = 0.1;
    
    
    [ p, q, ~ ] = size( I1 );
    ON = ones( p, q );

    % %structure check
    structureness = computeImageTexturedness( I1 );
    avgStructureness = mean( structureness( : ) );
    Off1 = structureness < max( structurePortion * avgStructureness, minValue );
    % Off1=zeros(size(ON));

    % Forward-backward check
    [ Ys_ori, Xs_ori ] = ndgrid( [ 1 : p ], [ 1 : q ] );
    Xs_new = Xs_ori + Ff( : , : , 1 );
    Ys_new = Ys_ori + Ff( : , : , 2 );
    flo_back_interp( : , : , 1 ) = interp2( Fb( : , : , 1 ), Xs_new, Ys_new );
    flo_back_interp( : , : , 2 ) = interp2( Fb( : , : , 2 ), Xs_new, Ys_new );
    Off2 = sqrt( ( Ff( : , : , 1 ) + flo_back_interp( : , : , 1 ) ) .^ 2 + ...
        ( Ff( : , : , 2 ) + flo_back_interp( : , : , 2 ) ) .^ 2 ) > ...
        fb_check * sqrt( Ff( : , : , 1 ) .^ 2 + Ff( : , : , 2 ) .^ 2 + ...
        flo_back_interp( : , : , 1 ) .^ 2 + ...
        flo_back_interp( : , : , 2 ) .^ 2 ) + fb_check2;
    
    if( debug > 3 )
        figure( 2 ),
        subplot( 2, 2, 1 );
        imshow( I1 )
        subplot( 2, 2, 2 );
        imagesc( sqrt( ( Ff( : , : , 1 ) + flo_back_interp( : , : , 1 ) ) .^ 2 + ...
            ( Ff( : , : , 2 ) + flo_back_interp( : , : , 2 ) ) .^ 2 ) );
        axis image;
        subplot( 2, 2, 3 );
        imagesc( sqrt( Ff( : , : , 1 ) .^ 2 + Ff( : , : , 2 ) .^ 2 + ...
            flo_back_interp( : , : , 1 ) .^ 2 + flo_back_interp( : , : , 2 ) .^ 2 ) );
        axis image;
        subplot( 2, 2, 4 );
        imagesc( sqrt( ( Ff( : , : , 1 ) + flo_back_interp( : , : , 1 ) ) .^2 + ...
            ( Ff( : , : , 2 ) + flo_back_interp( : , : , 2 ) ) .^ 2 ) > ...
            fb_check * sqrt( Ff( : , : , 1 ) .^ 2 + Ff( : , : , 2 ) .^ 2 + ...
            flo_back_interp( : , : , 1 ) .^ 2 + flo_back_interp( : , :, 2 ) .^ 2 ) + fb_check2 );
        axis image;
        title( 'forward backward check->high values indicate problem' )
    end

    ON = ON .* ( Off2 == 0 ) .* ( Off1 == 0 );

end

function structureness = computeImageTexturedness( im )

    if( size( im, 3 ) > 1 )
        im = im2double( rgb2gray( im ) );
    else
        im = im2double( im );
    end
    
    imageSmoothSigma = 3;
    gauss_1D = mkGaussian( imageSmoothSigma );
    gauss_2D = conv2( gauss_1D, gauss_1D' );
    dx = [ -1 0 1 ];   % Derivative masks
    dy = dx';
    gradient_filter_x = conv2( gauss_2D, dx, 'same' );
    gradient_filter_y = conv2( gauss_2D, dy, 'same' );
    % im_smooth=conv2(im,gauss_2D,'same');
    % Compute gradient
    Ix = conv2( im, gradient_filter_x, 'same' );      % Image derivatives
    Iy = conv2( im, gradient_filter_y, 'same' );
    imageSmoothSigma = 25;
    gauss_1D = mkGaussian( imageSmoothSigma );
    gauss_2D = conv2( gauss_1D, gauss_1D' );
    Ix2 = conv2( Ix .^ 2, gauss_2D, 'same' );
    Iy2 = conv2( Iy .^ 2, gauss_2D, 'same' );
    Ixy = conv2( Ix .* Iy, gauss_2D, 'same' );
    % Estimate the cornerness
    structureness = ( Ix2 .* Iy2 - Ixy .^ 2 ) ./ ( Ix2 + Iy2 + eps );

end

function gauss = mkGaussian( hsamples )
% creating 1D Gaussian, of size hsamples * 2 + 1
%
% the gaussian is truncated at x = +- tail, and there are samples samples
% inbetween, where samples = hsamples * 2 + 1
%
%  Jianbo Shi, 2007

    tail = 2;
    samples = hsamples * 2 + 1;

    x = linspace( - tail, tail, samples );
    gauss = exp( - x .^ 2 );
    %s = sum(gauss)/length(x);gauss = gauss-s;
    gauss = gauss / sum( abs( gauss ) );

    n = gauss * ones( samples, 1 );
    gauss = gauss / n;

end







