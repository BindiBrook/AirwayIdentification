clearvars; close all;clc

%% Open known airway files
% make list of file names
awpn='Data/ForAnalogRangeSelection/';
awfile={...
    'D34_1.1_13_Airways.mat',...
    'D34_1.2_13_Airways.mat',... % 'D34_1.3_13_Airways.mat',...
    'D34_1.4_13_Airways.mat',...
    'D34_2.1_13_Airways.mat',...
    'D34_2.2_13_Airways.mat',...
    'D34_4.1_13_Airways.mat',...
    'D34_4.2_13_Airways.mat',...
    'D34_4.3_13_Airways.mat',...
    'D34_4.4_13_Airways.mat',...
    'D34_5.1_13_Airways.mat',...
    'D34_5.2_13_Airways.mat',...
    'D35_5.3_13_Airways.mat',...
    'D35_6.1_13_Airways.mat',...
    'D35_6.2_13_Airways.mat',...
    'D35_6.3_13_Airways.mat',...
    'D35_6.4_13_Airways.mat',...
    'D35_7.1_13_Airways.mat',...
    };

% initialise variables
awfn_known=[]; awdensity_known=[]; selawtrans_known=[]; selawtheta_known=[];

another=1;
for another=1:length(awfile) % while another==1 % 
    clear awfn selawdensity selawtrans selawtheta % awfile awpn 
    % [awfile,awpn]=uigetfile('*.mat','Open file with known airways',awpn); 
        % get file location
    
    load([awpn awfile{another}]); % awfile]); % 
    curawfn=awfile{another}; % awfile; % 
    
    clear awfn_temp
    for kk=1:length(awfn) % Add prefix to filenames if missing
        awfn_temp=awfn{kk};
        if awfn_temp(1)=='s'
            awfn{kk}=[curawfn(1:10),'_',awfn_temp];
        end
    end
    
    awfn_known=[awfn_known awfn]; % get known airway file names
    awdensity_known=[awdensity_known selawdensity]; % get known airway pixel density in wall region
    selawtrans_known=[selawtrans_known selawtrans]; % get known airway number of transitions
    selawtheta_known=[selawtheta_known selawtheta]; % get known airway spoke angle differences
    
    % another=menu('Open another aiway file','yes','no');
end

%% Open all object files
% make list of file names
objpn='Data/ForAnalogRangeSelection/';
objfile={...
    'D34_1.1_13_Objects.mat',...
    'D34_1.2_13_Objects.mat',... % 'D34_1.3_13_Objects.mat',...
    'D34_1.4_13_Objects.mat',...
    'D34_2.1_13_Objects.mat',...
    'D34_2.2_13_Objects.mat',...
    'D34_4.1_13_Objects.mat',...
    'D34_4.2_13_Objects.mat',...
    'D34_4.3_13_Objects.mat',...
    'D34_4.4_13_Objects.mat',...
    'D34_5.1_13_Objects.mat',... 
    'D34_5.2_13_Objects.mat',...
    'D35_5.3_13_Objects.mat',...
    'D35_6.1_13_Objects.mat',...
    'D35_6.2_13_Objects.mat',...
    'D35_6.3_13_Objects.mat',...
    'D35_6.4_13_Objects.mat',...
    'D35_7.1_13_Objects.mat',...
    };
    
% initialise variables
objfn=[]; objdensity=[]; objtrans=[]; objtheta=[]; objDnum=[]; objMnum=[];

another=1;
for another=1:length(objfile) % while another==1 % 
    clear awfn selawdensity selawtrans selawtheta % objfile objpn 
    % [objfile,objpn]=uigetfile('*.mat','Open file with objects',objpn); 
        % get file location
    
    load([objpn objfile{another}]); % objfile]); % 
    curobjfn=objfile{another}; % objfile; % 
    
    objfn=[objfn awfn]; % get known airway file names
    objdensity=[objdensity selawdensity]; % get known airway pixel density in wall region
    objtrans=[objtrans selawtrans]; % get known airway number of transitions
    objtheta=[objtheta selawtheta]; % get known airway spoke angle differences
    
    % get "D" numbers and "M" or "mouse" numbers
    clear objDnum_var objMnum_var
    objDnum_var(1:length(awfn))=str2double(curobjfn(2:3));
    objMnum_var(1:length(awfn))=str2double(curobjfn(5:7));
    
    objDnum=[objDnum objDnum_var];
    objMnum=[objMnum objMnum_var];
    
    % another=menu('Open another object file','yes','no');
end

%% display means and ranges from Airways
disp(['selawtrans_known mean = ',num2str(mean(selawtrans_known)),' +/- std ',...
    num2str(std(selawtrans_known)), ' on range [',...
    num2str(min(selawtrans_known)), ',', ...
    num2str(max(selawtrans_known)),']'])

disp(['selawtheta_known mean = ',num2str(mean(selawtheta_known)),' +/- std ',...
    num2str(std(selawtheta_known)), ' on range [',...
    num2str(min(selawtheta_known)), ',', num2str(max(selawtheta_known)),']'])

disp(['awdensity_known mean = ',num2str(mean(awdensity_known)),' +/- std ',...
    num2str(std(awdensity_known)), ' on range [',...
    num2str(min(awdensity_known)), ',', num2str(max(awdensity_known)),']'])

%% display means and ranges from all Objects
disp(['objtrans mean = ',num2str(mean(objtrans)),' +/- std ',...
    num2str(std(objtrans)), ' on range [',...
    num2str(min(objtrans)), ',', ...
    num2str(max(objtrans)),']'])

disp(['objtheta mean = ',num2str(mean(objtheta)),' +/- std ',...
    num2str(std(objtheta)), ' on range [',...
    num2str(min(objtheta)), ',', num2str(max(objtheta)),']'])

disp(['objdensity mean = ',num2str(mean(objdensity)),' +/- std ',...
    num2str(std(objdensity)), ' on range [',...
    num2str(min(objdensity)), ',', num2str(max(objdensity)),']'])

%% Set ranges
trans_low=input('Enter lower limit on selawtrans ');
trans_up=input('Enter upper limit on selawtrans ');
theta_low=input('Enter lower limit on selawtheta ');
theta_up=input('Enter upper limit on selawtheta ');
dens_low=input('Enter lower limit on awdensity ');
dens_up=input('Enter upper limit on awdensity ');

%% First, check ranges of density and give score
% Get segment and object number from airway files
awDnum=zeros(1,length(awfn_known)); awMnum=zeros(1,length(awfn_known));
for aw=1:length(awfn_known) % iterate through known airway files
    awfn_str=awfn_known{aw}; % get current aiway file name
    % First, get segment number
    awDnum(aw)=str2double(awfn_str(2:3));
    awMnum(aw)=str2double(awfn_str(5:7));
    itn=1;
    clear aw_char
    for kk=20:length(awfn_str)
        if awfn_str(kk)=='_'
            break
        end
        aw_char(itn)=awfn_str(kk);
        itn=itn+1;
    end
    aw_seg(aw)=str2double(aw_char);
    % Then get Object number
    if aw_seg(aw)<10
        strt=35;
    else
        strt=36;
    end
    clear itn kk aw_char
    itn=1;
    for kk=strt:length(awfn_str)
        % disp(awfn_str(kk));
        if awfn_str(kk)=='.'
            break
        end
        aw_char(itn)=awfn_str(kk);
        itn=itn+1;
    end
    aw_num(aw)=str2double(aw_char);
end

% Get segment and object number from object files
for obj=1:length(objfn)
    objfn_str=objfn{obj};
    % First, get segment number
    itn=1;
    clear obj_char
    for kk=9:length(objfn_str)
        if objfn_str(kk)=='_'
            break
        end
        obj_char(itn)=objfn_str(kk);
        itn=itn+1;
    end
    obj_seg(obj)=str2double(obj_char);
    % Then get Object number
    if obj_seg(obj)<10
        strt=24;
    else
        strt=25;
    end
    clear itn kk obj_char
    itn=1;
    for kk=strt:length(objfn_str)
        if objfn_str(kk)=='.'
            break
        end
        obj_char(itn)=objfn_str(kk);
        itn=itn+1;
    end
    obj_num(obj)=str2double(obj_char);
end

% Iterate through objects and check if density is within range
dens_found=zeros(1,aw);
clear obj itn
itn=0; totfound=0;
for obj=1:length(objfn)
    % Check that density is within range
    if objdensity(obj)>=dens_low && objdensity(obj)<=dens_up && ...
        objtrans(obj)>=trans_low && objtrans(obj)<=trans_up && ...
        objtheta(obj)>=theta_low && objtheta(obj)<=theta_up
        totfound=totfound+1;
        % Now find the corresponding found airway
        clear kk
        for kk=1:aw
            if obj_seg(obj)==aw_seg(kk) && obj_num(obj)==aw_num(kk) && ...
                    objDnum(obj)==awDnum(kk) && objMnum(obj)==awMnum(kk)
                dens_found(kk)=1; % if matches, inter "1" into array
                itn=itn+1; % iterate to make array of matching file names
                awmatch{itn}=awfn_known{kk};
            end
        end
    end
end

% Display results
% disp(['Found a total of ',num2str(totfound),' objects']);
if sum(dens_found)==aw
    disp(['Found ',num2str(num2str(totfound)),' objects',...
        ' and all ',num2str(sum(dens_found)),' airways, or ',...
        num2str(100*sum(dens_found)/totfound),'% of searchable objects'])
else
    disp(['Found ',num2str(num2str(totfound)),' objects and ',...
        num2str(sum(dens_found)) ' of ',...
        num2str(aw),' airways, or ',...
        num2str(100*sum(dens_found)/totfound),'% of searchable objects'])
end
disp(['   compared to ',num2str(100*aw/obj),'% manually without this ',...
    'filtering, or ',num2str(aw),' airways out of ',num2str(obj),...
    ' total objects'])

%% Double check that all the airways are included with the objects
% Iterate through objects and check if density is within range
tot_found=zeros(1,aw);
clear obj itn
itn=0;
for obj=1:length(objfn)
    % Find the corresponding found airway
    clear kk
    for kk=1:aw
        if obj_seg(obj)==aw_seg(kk) && obj_num(obj)==aw_num(kk) && ...
                objDnum(obj)==awDnum(kk) && objMnum(obj)==awMnum(kk)
            tot_found(kk)=1; % if matches, inter "1" into array
            itn=itn+1; % iterate to make array of matching file names
            awtot{itn}=awfn_known{kk};
        end
    end
end
disp(['Check: found ',num2str(sum(tot_found)),' of ',num2str(aw),...
    ' airways in the object files'])