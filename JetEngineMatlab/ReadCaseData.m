function C = ReadCaseData(filename)
%READCASEDATA  Parse a group case-data file (e.g. Groep052.txt) into a struct.
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Expected file layout (one value per line, colon-separated):
%       Groep  52
%       Data
%                  Fuel :        H2
%              Tamb [K] :    300.00
%          P3overP2 [-] :      7.00
%             Pamb [Pa] : 100000.00
%        mfurate [kg/s] :      0.58
%                AF [-] :    170.35
%              v1 [m/s] :    200.00
%   Lines are matched by their label (case/whitespace/unit-insensitive),
%   so the exact column alignment and any [unit] tags don't matter.
%
%   Input:
%     filename - path to the case-data .txt file
%
%   Output:
%     C - struct with fields (all SI units, as given in the file):
%           Fuel      - fuel name string (must match an entry in Sp.Name)
%           Tamb      - ambient temperature [K]
%           P3overP2  - compressor pressure ratio [-]
%           Pamb      - ambient pressure [Pa]
%           mfurate   - fuel mass flow rate [kg/s]
%           AF        - air-to-fuel mass ratio [-]
%           v1        - inlet (flight) velocity [m/s]

fid = fopen(filename,'r');
if fid == -1
    error('ReadCaseData:FileNotFound','Could not open case file: %s',filename);
end
raw = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines = raw{1};

C = struct('Fuel',[],'Tamb',[],'P3overP2',[],'Pamb',[],'mfurate',[],'AF',[],'v1',[]);
for k = 1:numel(lines)
    line = lines{k};
    if ~contains(line,':')
        continue
    end
    idx   = find(line==':',1,'first');
    key   = line(1:idx-1);
    value = strtrim(line(idx+1:end));

    key = regexprep(key,'\[[^\]]*\]','');                                  % drop any [unit] tag
    key = lower(regexprep(strtrim(key),'\s+',''));                         % normalise: lowercase, no spaces

    switch key
        case 'fuel',      C.Fuel     = value;
        case 'tamb',      C.Tamb     = str2double(value);
        case 'p3overp2',  C.P3overP2 = str2double(value);
        case 'pamb',      C.Pamb     = str2double(value);
        case 'mfurate',   C.mfurate  = str2double(value);
        case 'af',        C.AF       = str2double(value);
        case 'v1',        C.v1       = str2double(value);
    end
end

reqFields = {'Fuel','Tamb','P3overP2','Pamb','mfurate','AF','v1'};
for i=1:numel(reqFields)
    if isempty(C.(reqFields{i}))
        error('ReadCaseData:MissingField','Field "%s" was not found in %s',reqFields{i},filename);
    end
end
end
