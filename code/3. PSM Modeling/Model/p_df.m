function [ ret ] = p_df( no )
global SYMBOL_VERIFICATION;

    if SYMBOL_VERIFICATION
        eval(['syms p_d', num2str(no), ';']);
        if no~=1
            ret = eval(['p_d', num2str(no),';']);
        else
            ret = 0;
        end
    else
        ret = p_bf(no)+p_cf(no);
    end
end



