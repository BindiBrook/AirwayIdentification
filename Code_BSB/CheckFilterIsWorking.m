clearvars; close all; clc

%% First, get info from KNOWN airway file names
% Set paths for known airways
knownawpn='/local/pmzmrh/r-drive/Mike/HistomorphologyMATLABCode/OVA13/';
specimen={...
    'D34_1.1_13',... % 1
    'D34_1.2_13',... % 2
    'D34_1.3_13',... % 3
    'D34_1.4_13',... % 4
    'D34_2.1_13',... % 5
    'D34_2.2_13',... % 6
    'D34_4.1_13',... % 7
    'D34_4.2_13',... % 8
    'D34_4.3_13',... % 9
    'D34_4.4_13',... % 10
    'D34_5.1_13',... % 11
    'D34_5.2_13',... % 12
    'D35_5.3_13',... % 13
    'D35_6.1_13',... % 14
    'D35_6.2_13',... % 15
    'D35_6.2_14',... % 16
    'D35_6.3_13',... % 17
    'D35_6.4_13',... % 18
    'D35_7.1_13'};   % 19

% Set path for directory containing objects for comparison
filterpn='/local/pmzmrh/';
folder1pn='Folder1';
folder2pn='Folder2';

diary([knownawpn 'GroundTruthTest.log'])
diary on

% Loop through specimens
for jj=1:length(specimen)
    disp(['Data for ' specimen{jj}])
    
    clear subpn
    subpn=[specimen{jj} '/Airways/'];
    
    % Get segment and object number from airway files
    clear awfolderfiles nawfiles
    awfolderfiles = dir([knownawpn subpn]); % all files in folder
    nawfiles = length(awfolderfiles);    % Number of files in folder
    aw=0;
    for ii=1:nawfiles % iterate through known airway files
        clear currentfilename
        currentfilename=awfolderfiles(ii).name;
        if length(currentfilename)>4 && currentfilename(1)~='.'
            if currentfilename(length(currentfilename)-3:...
                    length(currentfilename))=='.mat'
                finalword=currentfilename(length(currentfilename)-8:...
                    length(currentfilename)-4);
                if finalword(1)~='f' && finalword(2)~='i' && ...
                        finalword(3)~='n' && finalword(4)~='a' && ...
                        finalword(5)~='l'
                    aw=aw+1;
                    awfn_str=currentfilename; % get current aiway file name
                    if currentfilename(1)=='D'
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
                    else
                        % First, get segment number
                        % Get object "D" and "M" numbers
                        thisspecimen=specimen{jj};
                        awDnum(aw)=str2double(thisspecimen(2:3));
                        awMnum(aw)=str2double(thisspecimen(5:7));
                        itn=1;
                        clear aw_char
                        for kk=9:length(awfn_str)
                            if awfn_str(kk)=='_'
                                break
                            end
                            aw_char(itn)=awfn_str(kk);
                            itn=itn+1;
                        end
                        aw_seg(aw)=str2double(aw_char);
                        % Then get Object number
                        if aw_seg(aw)<10
                            strt=24;
                        else
                            strt=25;
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
                end
            end
        end
    end
    
    % Next, get info from Objects sorted by Filter
    % Set paths for filtered objects
    clear sfubpn folder1files folder2files nf1files nf2files
    if jj==1 || jj==7
        fsubpn=[specimen{jj} '/objects/'];
    else
        fsubpn=[specimen{jj} '/'];
    end
    
    % Get segment and object number from object files
    clear folder1files folder2files nf1files nf2files
    folder1files = dir([filterpn fsubpn folder1pn]); % all files in folder
    folder2files = dir([filterpn fsubpn folder2pn]); % all files in folder
    nf1files = length(folder1files); % Number of files in folder
    nf2files = length(folder2files); % Number of files in folder
    
    obj=0; % initialize object counter
    
    % FOLDER 1
    % Get segment and object number from object files
    for ii=1:nf1files
        clear objfn_str
        objfn_str=folder1files(ii).name;
        if length(objfn_str)>4 && objfn_str(1)~='.'
            obj=obj+1;
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
            
            % Get object "D" and "M" numbers
            objDnum(obj)=str2num(fsubpn(2:3));
            objMnum(obj)=str2num(fsubpn(5:7));
        end
    end
    
    % FOLDER 2
    obj2=0; % initilise counter
    % Get segment and object number from object files
    for ii=1:nf2files
        clear objfn_str
        objfn_str=folder2files(ii).name;
        if length(objfn_str)>4 && objfn_str(1)~='.'
            obj2=obj2+1;
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
            obj2_seg(obj2)=str2double(obj_char);
            % Then get Object number
            if obj2_seg(obj2)<10
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
            obj2_num(obj2)=str2double(obj_char);
            
            % Get object "D" and "M" numbers
            obj2Dnum(obj2)=str2num(fsubpn(2:3));
            obj2Mnum(obj2)=str2num(fsubpn(5:7));
        end
    end
    
    % Iterate through objects and check if airways are in this set of objects
    tot_found=zeros(1,aw);
    n=obj; clear obj itn
    itn=0;
    for obj=1:n
        % Find the corresponding found airway
        clear kk
        for kk=1:aw
            if obj_seg(obj)==aw_seg(kk) && obj_num(obj)==aw_num(kk) && ...
                    objDnum(obj)==awDnum(kk) && objMnum(obj)==awMnum(kk)
                tot_found(kk)=1; % if matches, inter "1" into array
                itn=itn+1; % iterate to make array of matching file names
                % awtot{itn}=awfn_known{kk};
            end
        end
    end
    disp(['Check: found ',num2str(sum(tot_found)),' of ',num2str(aw),...
        ' airways in object files in Folder 1'])
    
    % Iterate through objects and check if airways are in this set of objects
    tot_found=zeros(1,aw);
    clear n; n=obj2; clear obj itn
    itn=0;
    for obj=1:n
        % Find the corresponding found airway
        clear kk
        for kk=1:aw
            if obj2_seg(obj)==aw_seg(kk) && obj2_num(obj)==aw_num(kk) && ...
                    obj2Dnum(obj)==awDnum(kk) && obj2Mnum(obj)==awMnum(kk)
                tot_found(kk)=1; % if matches, inter "1" into array
                itn=itn+1; % iterate to make array of matching file names
                % awtot{itn}=awfn_known{kk};
            end
        end
    end
    disp(['Check: found ',num2str(sum(tot_found)),' of ',num2str(aw),...
        ' airways in object files in Folder 2'])
    
    disp('--'); % make a space
    
end
diary off