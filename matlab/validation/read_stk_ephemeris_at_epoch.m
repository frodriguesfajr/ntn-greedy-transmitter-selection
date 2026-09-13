function out = read_stk_ephemeris_at_epoch(ephFile,epochUTC)
%READ_STK_EPHEMERIS_AT_EPOCH Read STK Fixed-frame state at an exact epoch.

if ~isfile(ephFile)
    error('STK ephemeris not found: %s',ephFile);
end

fid = fopen(ephFile,'r');
if fid < 0
    error('Could not open: %s',ephFile);
end

cleanup = onCleanup(@() fclose(fid));

scenarioEpoch = [];
dataStart = false;

while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end

    s = strtrim(line);

    if startsWith(s,'ScenarioEpoch ')
        epochText = strtrim(extractAfter(string(s),'ScenarioEpoch '));
        scenarioEpoch = datetime(epochText, ...
            'InputFormat','dd MMM yyyy HH:mm:ss.SSS', ...
            'Locale','en_US','TimeZone','UTC');
    end

    if strcmp(s,'EphemerisTimePosVel')
        dataStart = true;
        break
    end
end

if isempty(scenarioEpoch)
    error('ScenarioEpoch not found in %s',ephFile);
end

if ~dataStart
    error('EphemerisTimePosVel not found in %s',ephFile);
end

C = textscan(fid,'%f %f %f %f %f %f %f');
data = [C{:}];

if isempty(data)
    error('No ephemeris samples found in %s',ephFile);
end

if isempty(epochUTC.TimeZone)
    epochUTC.TimeZone = 'UTC';
end

target_s = seconds(epochUTC - scenarioEpoch);

[sampleError_s,idx] = min(abs(data(:,1)-target_s));

if sampleError_s > 1e-6
    error('Requested epoch is not an exact STK sample. Nearest offset = %.6f s.',sampleError_s);
end

out.scenarioEpoch = scenarioEpoch;
out.offset_s       = data(idx,1);
out.rECEF_m        = data(idx,2:4).';
out.vECEF_mps      = data(idx,5:7).';

end
