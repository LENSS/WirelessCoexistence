function [ ret ] = p_bf( no )
global SYMBOL_VERIFICATION;

    if SYMBOL_VERIFICATION

        eval(['syms p_b', num2str(no), ';']);
        if no~=1
            ret = eval(['p_b', num2str(no),';']);
        else
            ret = 0;
        end
    
    else
        ret = 1-p_af(no)-p_cf(no);
        
    end
end

