clear;

%runDate = '2021-01-07';
%folderToImport = ['/Users/bramsona/Dropbox/Science/ThermalModel/NPLD_Lag_CollabWithAlyssa/v12/Profile2/Wind/Wind20/Output/', runDate, '/'];
%fileToImport = dir([folderToImport, '*_output.txt']);
fileToImport = dir([pwd, '/', '*_output.txt']);
stname = fileToImport(1).name;
output = importfile([pwd, '/', stname]);

%assignin('base', ['output_', slopeStr], output); % Change name
%clear output;

FFsub_tab = output(:,10);
LTime_tab = output(:,1);

FFsubl = table2array(FFsub_tab(:,:));
LaskarTimestep = table2array(LTime_tab(:,:));
% Not using h0s right now so set them all to -1
h0s = FFsubl*0 -1;

save(['Layer2W_500ka.mat'], 'FFsubl', 'h0s', 'LaskarTimestep');