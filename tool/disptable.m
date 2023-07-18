function disptable(colname,varargin)
%% DISPTABLE display data in table
%
% DISPTABLE({'Data 1', 'Data 2'},data1,format1,data2,format2) display data1
% and data2 in table using format1 and format2, respectively.  
%
% Example
%
%	load dispdata
%   colname = {'#nodes'     '|u_I-u_h|_1'      '||p-p_h||'};
%   disptable(colname,err.N,[],err.uH1,'%0.5e',err.pL2,'%0.5e');
%
% See also displaytable
%
% Copyright (C)  Long Chen. See COPYRIGHT.txt for details.

% displaytable(colname,varargin{:});
col = varargin{1};
rownum = size(col,1);
colnum = size(colname,2);
colwidth = zeros(1,colnum);
dispresult = char(zeros(rownum,6*colnum));
idx = 1;
for c = 1:nargin/2
    if isempty(varargin{2*c})
        currentstr = num2str(varargin{2*c-1});
    else
        currentstr = num2str(varargin{2*c-1},varargin{2*c});
    end
    colwidth(c) = size(currentstr,2);
    dispresult(:,idx:idx+colwidth(c)-1) = currentstr; 
    idx = idx+colwidth(c)+3;
end
for c = 1:colnum
    spnum = colwidth(c) - length(colname{c})+3;
    if spnum >0
        if c>1
            colname{c} = [char(zeros(1,floor(spnum/2))) colname{c} ...
                          char(zeros(1,ceil(spnum/2)))];
        else
            colname{c} = [char(zeros(1,floor(spnum/3))) colname{c} ...
                          char(zeros(1,ceil(spnum/2)))];            
        end
    end
end
% fprintf('\n');
disp(cell2mat(colname));
fprintf('\n');
disp(dispresult);
fprintf('\n');