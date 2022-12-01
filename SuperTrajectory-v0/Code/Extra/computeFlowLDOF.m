function [fflow, bflow] = computeFlowLDOF(data, options)
    flowfolder = fullfile(options.datafolder, 'flow');
    if( ~exist( flowfolder, 'dir' ) ) 
         mkdir( flowfolder ); 
    end;
    
    for framenum = 1:data.nframe-1
        fprintf( 'processing: Frame #%i\n', framenum );
        filename = fullfile( flowfolder, ['foreflow_' data.names{framenum} '.mat'] );
        if(exist( filename, 'file' ) )
            load( filename);          
        else       
            foreflow = mex_LDOF( double(data.frames{framenum}), double(data.frames{framenum+1}) );
            save( filename, 'foreflow', '-v7.3' ); 
        end
        fflow{framenum} = foreflow;
        filename = fullfile( flowfolder, ['backflow_' data.names{framenum} '.mat'] );
        if(exist( filename, 'file' ) )
            load( filename);          
        else 
            backflow = mex_LDOF( double(data.frames{framenum+1}), double(data.frames{framenum}) );
            file = fullfile( flowfolder, ['backflow_' data.names{framenum} '.mat'] );
            save( filename, 'backflow', '-v7.3' ); 
        end
        bflow{framenum} = backflow;
     end
end
