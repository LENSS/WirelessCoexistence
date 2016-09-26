function [ ret ] = p_af( no )
global SYMBOL_VERIFICATION;
global tauc;

    if SYMBOL_VERIFICATION

        eval(['syms p_a', num2str(no), ';']);
        if no~=1
            ret = eval(['p_a', num2str(no),';']);
        else
            ret = 1;
        end
    else
        
        ret = tauc^(no-1);
        
    end
end

%%%%%