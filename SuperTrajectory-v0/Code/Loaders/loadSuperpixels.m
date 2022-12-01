function [superpixelsLAB,superpixelsLabel]= loadSuperpixels( options )
    superpixelsLAB = [];
    superpixelsLabel = [];
    file = fullfile( options.datafolder, 'superpixels.mat' );
    if( exist( file, 'file' ) )
        superpixel = load( file ); 
        superpixelsLAB = superpixel.superpixelsLAB;           
        superpixelsLabel = superpixel.superpixelsLabel;
    end
    
end