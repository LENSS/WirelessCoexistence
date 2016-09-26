function [ out ] = zw_printf( varargin )
global DEBUG;
%     fprintf('Total number of inputs = %d\n', nargin);
%     celldisp(varargin);
    if DEBUG
        out = fprintf(varargin{:});
    else
        out = 0;
    end
end

% function mfprintf(fid, varargin)
%     arrayfun(@(fid) fprintf(fid, varargin{:}), fid);
% end