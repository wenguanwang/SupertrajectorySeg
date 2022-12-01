
function [ XYTRGB, tr_id ] = quick_tr( tr )

    if isempty( tr ) 
        XYTRGB = []; 
        tr_id = []; 
        tr = []; 
        return; 
    end
    
    XYTRGB = cat( 2, tr.XYTRGB );
    tr_id = zeros( 1, size( XYTRGB, 2 ) );
    c = 1;
    for ii = 1 : length( tr )
        tr_id( c : c + size( tr( ii ).XYTRGB, 2 ) - 1 ) = ii;
        c = c + size( tr( ii ).XYTRGB, 2 );
    end
end
