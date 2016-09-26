function [ out ] = zw_printf( varargin )
global DEBUG_PRT DEBUG_LOG;
%     fprintf('Total number of inputs = %d\n', nargin);
%     celldisp(varargin);
    if DEBUG_PRT
        
        if DEBUG_LOG; diary on; end
        out = fprintf(varargin{:});
        if DEBUG_LOG; diary off; end
    else
        out = 0;
    end
end

% function mfprintf(fid, varargin)
%     arrayfun(@(fid) fprintf(fid, varargin{:}), fid);
% end