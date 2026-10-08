function PrintVar(name,value,unit)
%PRINTVAR  Print one quantity on its own line: name, value, unit (aligned)
%   Used by Assignment.m and all stages so the console output is a vertical
%   list with one value per line. value may also be text (e.g. the method).
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Input:
%     name   - label to print
%     value  - number or text
%     unit   - unit string (optional)

if nargin < 3, unit = ''; end
if ischar(value) || isstring(value)
    fprintf('  %-22s %14s\n',name,value);
else
    fprintf('  %-22s %14.4f  %s\n',name,value,unit);
end
end
