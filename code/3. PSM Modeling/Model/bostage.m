function [ ret ] = bostage( index )
%return the backoff stage
%the smallest stage is 1 ,not 0

global cwarray;
global snarray;

    if index<=0
        ret = index;
        return;
    end

    ret=find(snarray>=index,1);
    if isempty(ret)
        
        ret=-2;
        return;
        
    end    
%     disp(idx);
%     if idx~=1
%         ret = index-cwarray(idx-1);
%     else
%         ret = index;
%     end

end

