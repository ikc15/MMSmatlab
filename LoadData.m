function [fgm, edp, fpii, fpie]=LoadData(instrument)
global yyyymmddHHMM starttime endtime cavitystart cavityend; 
%% INPUTS
yyyymmddHHMM = '201601210106'; % enter date and time e.g. 19th Jan 2017 at 00:33 is 201701190033
starttime = '01:06:41.10'; % start of current sheet
endtime = '01:06:52.04'; % end of current sheet
cavitystart = '01:06:50.40';
cavityend = '01:06:51.90';
cavity = 1;    % 0 for no cavity
if cavity == 0    
    cavitystart = starttime;
    cavityend = endtime;
end
folder = '/Users/Andrea/Documents/MSci Project/Data'; % location of data folder

%% 
files = dir(fullfile(folder,sprintf('**/*%s*',yyyymmddHHMM)));

if ismember('fgm', instrument)
    %Flux gate magnetometer
    fgm = spdfcdfread(strcat(files(2).folder,'/',files(2).name),'ShowProgress',true);%, 'KeepEpochAsIs',true);
end
if ismember('edp', instrument)
    % electric field double probe
    edp = spdfcdfread(strcat(files(1).folder,'/',files(1).name), 'ShowProgress', true);
end
if ismember('fpi', instrument)
    % Fast plasma instrument
    fpie = spdfcdfread(strcat(files(3).folder,'/',files(3).name), 'ShowProgress', true);
    fpii = spdfcdfread(strcat(files(4).folder,'/',files(4).name), 'ShowProgress', true);
end
end