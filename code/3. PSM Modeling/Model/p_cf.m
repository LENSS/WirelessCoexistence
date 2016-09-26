function [ ret ] = p_cf( no )
global SYMBOL_VERIFICATION;
global tau;
global tauc;

    if SYMBOL_VERIFICATION
        eval(['syms p_c', num2str(no), ';']);
        if no~=1
            ret = eval(['p_c', num2str(no),';']);
        else
            ret = 0;
        end
    else
        if no==1
            ret = 0;
        else
            ret = (no-1)*tau*tauc^(no-2);
        end
    end
end
%%%%%%%%
