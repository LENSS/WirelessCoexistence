function [ ret ] = compute_comb(zw)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global nn;
    ret = zeros(nn,nn); 
    cc=0;
    for n=nn+1-zw:-1:1
        ccc=0;
        for k=1:nn-cc
%                     zw_printf('%d, %d\n', n-1,cc+k-1);
%                     fprintf('%d, %d\n', nn-cc-1,ccc);
            ret(nn-cc,ccc+1)=nchoosek(nn-cc-1,ccc);
%                 zw_printf('%d to %d, %d*P^%d*(1-P)^%d\n', n-1, cc+k, nchoosek(nn-cc-1,ccc), nn-cc-1-ccc, ccc);
            ccc=ccc+1;
        end
        cc=cc+1;
    end
%                     fprintf('\n');

end

