clearvars; clc; close all;

warning('off','all'); % turn off warnings;

% [fn,pn]=uigetfile('Open image'); % get file location

mainpn='/local/pmzmrh/r-drive/Mike/HistomorphologyMATLABcode/OVA13/';

subpn_list={...
    'D35_5.3_13/Airways/',...
    'D35_5.3_13/Objects/',...
    'D35_6.1_13/Airways/',...
    'D35_6.1_13/Objects/',...
    'D35_6.2_13/Airways/',...
    'D35_6.2_13/Objects/',...
    'D35_6.2_14/Airways/',...
    'D35_6.2_14/Objects/',...
    'D35_6.3_13/Airways/',...
    'D35_6.3_13/Objects/',...
    'D35_6.4_13/Airways/',...
    'D35_6.4_13/Objects/',...
    'D35_7.1_13/Airways/',...
    'D35_7.1_13/Objects/',...
    };

%     'D34_1.1_13/Airways/',...
%     'D34_1.1_13/Non-Airway_Objects/',...
%     'D34_1.2_13/Airways/',...
%     'D34_1.2_13/Non-Airway_Objects/',...
%     'D34_1.3_13/Airways/',...
%     'D34_1.3_13/Non-Airway_Objects/',...
%     'D34_1.4_13/Airways/',...
%     'D34_1.4_13/Non-Airway_Objects/',...
%     'D34_2.1_13/Airways/',...
%     'D34_2.1_13/Non-Airway_Objects/',...
%     'D34_2.2_13/Airways/',...
%     'D34_2.2_13/Non-Airway_Objects/',...
%     'D34_4.1_13/Airways/',...
%     'D34_4.1_13/Objects/',...
%     'D34_4.2_13/Airways/',...
%     'D34_4.2_13/Objects/',...
%     'D34_4.3_13/Airways/',...
%     'D34_4.3_13/Objects/',...
%     'D34_4.4_13/Airways/',...
%     'D34_4.4_13/Objects/',...
%     'D34_5.1_13/Airways/',...
%     'D34_5.1_13/Objects/',...
%     'D34_5.2_13/Airways/',...
%     'D34_5.2_13/Objects/',...

for folderitn=1:length(subpn_list)
    % clear variables every loop
    clear subpn pn folderfiles nfiles itn awfnitn
    clear currentfilename finalword tottrans meanTheta awdensity flag
    clear awfn selawtrans selawtheta selawdensity
    
    subpn=subpn_list{folderitn};
    
    % 'D34_1.1_13/Non-Airway_Objects/'; % Airways/';
    % input('Enter folder and subfolder name '); %
    
    pn=[mainpn subpn];
    
    % Get all JPG files in this directory
    folderfiles = dir(pn); % all files in folder
    nfiles = length(folderfiles);    % Number of files found
    itn=1; awfnitn=0;
    for ii=1:nfiles
        currentfilename = folderfiles(ii).name;
        if length(currentfilename)>4 && currentfilename(1)~='.'
            if currentfilename(length(currentfilename)-3:...
                    length(currentfilename))=='.mat'
                finalword=currentfilename(length(currentfilename)-8:...
                    length(currentfilename)-4);
                if finalword(1)~='f' && finalword(2)~='i' && ...
                        finalword(3)~='n' && finalword(4)~='a' && ...
                        finalword(5)~='l'
                    [tottrans(itn),meanTheta(itn),awdensity(itn),flag]=...
                        TracePixelsAlongLine(pn,currentfilename);
                    if isnan(tottrans(itn))~=1 && ...
                            isnan(meanTheta(itn))~=1 && ...
                            isnan(awdensity(itn))~=1
                        awfnitn=awfnitn+1;
                        awfn{awfnitn}=currentfilename;
                        selawtrans(awfnitn)=tottrans(itn);
                        selawtheta(awfnitn)=meanTheta(itn);
                        selawdensity(awfnitn)=awdensity(itn);
                    end
                end
            end
        end
    end
    
    %% Save output to csv file
    % Clear variables in this loop
    clear svfn
    
    % Create, populate, and save table of values
    if subpn(12)=='N'
        svfn=[subpn(1:10) '_Objects'];
    elseif subpn(12)=='A' || subpn(12)=='O'
        svfn=[subpn(1:10) '_' subpn(12:18)];
    end
    svpnfn=[svfn,'.csv']; % make csv file name
    ColLabels={'filename',...
        'awdensity',...
        'tottrans',...
        'meanTheta'};
    
    if exist(svpnfn,'file')
        disp([svpnfn,' exists, will append data']);
    else
        disp([svpnfn,' does not exist, will create file']);
        fileID=fopen(svpnfn,'w');
        fprintf(fileID,'%s %s %s %s %s %s %s \n',...
            ColLabels{1},',',ColLabels{2},',',ColLabels{3},',',...
            ColLabels{4});
        fclose(fileID); clear fileID
    end
    
    fileID=fopen(svpnfn,'a');
    for kk=1:awfnitn
        fprintf(fileID,'%s %s %4.8f %s %4.8f %s %4.8f \n',...
            awfn{kk},',',...
            selawdensity(kk),',',...
            selawtrans(kk),',',...
            selawtheta(kk));
    end
    dispmess=msgbox('Process Complete');
    
    svpnfn=[svfn,'.mat']; % make mat file name
    save(svpnfn,'awfn','selawdensity','selawtrans','selawtheta');
end


% ind= ~isnan(tottrans); tottrans=tottrans(ind);
% disp(['totrans mean = ',num2str(mean(tottrans)),' +/- std ',...
%     num2str(std(tottrans)), ' on range [',...
%     num2str(min(tottrans)), ',', num2str(max(tottrans)),']'])
%
% ind= ~isnan(meanTheta); meanTheta=meanTheta(ind);
% disp(['meanTheta mean = ',num2str(mean(meanTheta)),' +/- std ',...
%     num2str(std(meanTheta)), ' on range [',...
%     num2str(min(meanTheta)), ',', num2str(max(meanTheta)),']'])
%
% ind= ~isnan(awdensity); awdensity=awdensity(ind);
% disp(['awdensity mean = ',num2str(mean(awdensity)),' +/- std ',...
%     num2str(std(awdensity)), ' on range [',...
%     num2str(min(awdensity)), ',', num2str(max(awdensity)),']'])

% % Check tottrans is within range
% translim=[0 100]; thetalim=[0 90]; densitylim=[75 150];
% if tottrans(itn)>translim(1) && ...
%         tottrans(itn)<translim(2)
%     awtrans(awfnitn)=1;
% else
%     awtrans(awfnitn)=0;
% end
% % Check meanTheta is within range
% if meanTheta(itn)>thetalim(1) && ...
%         meanTheta(itn)>thetalim(2)
%     awtheta(awfnitn)=1;
% else
%     awtheta(awfnitn)=0;
% end
% % Check awdensity is within range
% if awdensity(itn)>densitylim(1) && ...
%         awdensity(itn)>densitylim(2)
%     selawdens(awfnitn)=1;
% else
%     selawdens(awfnitn)=0;
% end