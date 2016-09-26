function [ ret ] = windowleft( index )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

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
    
    
    ret = mod(ret,(nn*windowsize));
    if ret==0
        ret=1;
        return;
    end
    t = floor(ret/nn);
    
    if mod(ret,nn)==0
        t=t-1;
    end
    ret = windowsize-t;


end

