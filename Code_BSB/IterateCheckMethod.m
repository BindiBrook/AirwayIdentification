clearvars; close all; clc

%% File Input
pn='/Volumes/R11/Maths-RespMedCollab/Mike/HistomorphologyMATLABCode/OVA13/';
Airways_Folder1=1; % Find duplicates in "Airways" or "Folder1"?

fnarray={...
    'D37_5.4_13',...
    };
    
%     'D41_20.3_1',...
%     'D41_18.2_1',...
%     'D41_18.3_1',...
%     'D41_18.4_1',...
%     'D41_20.1_1',...
%     'D41_20.2_1',...
%     'D41_20.4_1',...

for kk=1:length(fnarray)
    CheckForDuplicates([pn fnarray{kk} '/'], 'InfoFile.mat',Airways_Folder1)
end

% for folderitn=1:length(subpn_list)
%     % clear variables every loop
%     clear subpn pn folderfiles nfiles itn awfnitn
%     clear currentfilename finalword tottrans meanTheta awdensity flag
%     clear awfn selawtrans selawtheta selawdensity
%     
%     subpn=subpn_list{folderitn};