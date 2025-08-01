%% Header

filename = filenamelist(z).name;
header = struct;
fid = fopen([path_data,filename],'r','ieee-le');
fgetl(fid);
for i = 1:100
    sv = fgetl(fid); % Read the data in this row
    ind = isstrprop(sv,'digit'); % Find numbers in a string
    a = find(ind > 0,1);
    if isempty(a)
        ind2 = isstrprop(sv,'wspace');
        a = find(ind2>0,1);
        if isempty(a)
            a = length(sv)+1;
        end
    end
    sv2 = sv(1:a-1);
    inda = isstrprop(sv2,'alpha'); % Find numbers in a string
    svar = sv2(inda); % The data text
    eval(['header.',svar,'=char(sv(a:end));']);
    if strcmpi(svar,'data')
        break
    end
end
%
if isfield(header,'GPSSTATUSLOCK')
    ind_gps = find(isstrprop(header.GPSSTATUSLOCK,'alpha') > 0);
    lat = str2double(header.GPSSTATUSLOCK(1:ind_gps(1)-2));
    lon = str2double(header.GPSSTATUSLOCK(ind_gps(1)+2:ind_gps(2)-2));
else
    lat = 25.0;
    lon = 121.1;
end

time1 = str2double(header.DATA);
ut = floor(lon/15);
unix_time_start = datenum('1-jan-1970');                % days of 0000/00/00~1969/12/31
dateIndex = (time1+60*60*8)/60/60/24+unix_time_start;   % days of 0000/00/00~data time
[year0,mon0,day0,h0,m0,s0] = datevec(dateIndex);
Nchannel = str2double(header.CHANNELS);
freq_dwell = str2double(header.FREQUENCYDWELL);
freq_start = str2double(header.FREQUENCYSTART);
freq_step = str2double(header.FREQUENCYSTEP);
freq_stop = str2double(header.FREQUENCYSTOP);
freq = freq_start:freq_step:freq_stop;
% freq=freq(1:161);
pnum = str2double(header.PULSENUMBER);
ind_PULSELENGTH = find(isstrprop(header.PULSELENGTH,'wspace')>0);

if isempty(ind_PULSELENGTH),ind_PULSELENGTH = length(header.PULSELENGTH)+1;end

pw = str2double(header.PULSELENGTH(1:ind_PULSELENGTH-1))*1e-9;
ind_CODE = find(isstrprop(header.PULSECODE,'wspace')>0);

if isempty(ind_CODE),ind_CODE = length(header.PULSECODE)+1;end

codelength = length(header.PULSECODE(1:ind_CODE-1));
nrg1 = str2double(header.GATES);
cint = str2double(header.INTEGRATIONS);
min_range = str2double(header.RANGE)/1000;
rench = str2double(header.CHANNELS)/2;          % num of receiver
r_res = pw*3e8/2/1000;                          % range resolution (unit:km)
max_range = min_range+r_res*(nrg1-1);           % min_range and max_range (unit:km)
range = min_range+r_res:r_res:max_range;
ntrial = freq_dwell*length(freq)/cint;
fSampleLength = nrg1*rench*ntrial;
data = fread(fid, fSampleLength,'short');
fclose(fid);
year_date_time = (datestr(dateIndex,31));

if length(data)~=fSampleLength
    lim = length(data);
    data(lim+1:fSampleLength) = nan;
end