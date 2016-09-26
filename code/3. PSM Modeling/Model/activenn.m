function [ ret ] = activenn( index )
%return the number of active nodes from index
global nn;

    if index<=0
        ret = index;
        return;
    end
    
	ret=nn+1-mod(index,nn);
	if ret==nn+1;
		ret=1;
	end
end

