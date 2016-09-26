function [ pex ] = compute_pxn( zw, P )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    global nn;
    global comb_container;
    
    pex = zeros(nn,nn); 
    cc=0;
    for n=nn+1-zw:-1:1
        ccc=0;
        for k=1:nn-cc
%                     zw_printf('%d, %d\n', n-1,cc+k-1);
%                     zw_printf('%d, %d\n', nn-cc-1,ccc);
%             pex(n,cc+k)=nchoosek(nn-cc-1,ccc)*P^(nn-cc-1-ccc)*(1-P)^(ccc);
            pex(n,cc+k)=comb_container{zw}(nn-cc,ccc+1)*P^(nn-cc-1-ccc)*(1-P)^(ccc);
%                 zw_printf('%d to %d, %d*P^%d*(1-P)^%d\n', n-1, cc+k, nchoosek(nn-cc-1,ccc), nn-cc-1-ccc, ccc);
            ccc=ccc+1;
        end
        cc=cc+1;
    end
%     zw_printf('pex:\n');
%     zw_disp(pex);

end

