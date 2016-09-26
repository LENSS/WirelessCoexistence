function [ ret ] = bocounter( index )
%return the backoff counter, from cw to 1 (not cw-1 to 0)

global cwarray;
global snarray;
global windowsize;
global nn;

    if index<=0
        ret = index;
        return;
    end

    idx=find(snarray>=index,1);
    
    if isempty(idx)
        
        ret=-2;
        return;
        
    end
%     disp(idx);
    if idx~=1
        ret = index-snarray(idx-1);
    else
        ret = index;
    end

%     fprintf('1 ret: %d\n', ret);

    if mod(ret,(nn*windowsize))==0
        ret = floor(ret/(nn*windowsize))-1;
        ret=cwarray(idx)-ret;
%         fprintf('3 ret: %d\n', ret);
        return;
    end

    ret = floor(ret/(nn*windowsize));
%     disp(ret);
    ret = cwarray(idx)-ret;
%         fprintf('2 ret: %d\n', ret);
% ret = 0;
    end



