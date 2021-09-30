clearvars; close all; clc;
% Program to select which objects are actually airways


%% File Input
dispmess=msgbox(['Navigate to and select the tif file you were analysing.']);
drawnow; waitfor(dispmess); clear dispmess

[file,path]=uigetfile('../*.tif','Select tif file you were analysing');

% Check to see if there is a sub-directory "Airways" and create if not
% Check to see if the "log file" containing airway ID exists, and create if not
pnAirways = [path file(1:end-4) '/Airways'];
IDlogpath = [path file(1:end-4) '/IDlog.mat'];
if exist(pnAirways)
    if exist(IDlogpath)
        load(IDlogpath)
        startID=lastID+1;
    else %it's not clear that this is needed - if Airways directory exists, so will the log? But put in anyway.
        lastID = 2; startID=lastID+1;
        save(IDlogpath,'lastID');        
    end
else
    mkdir(pnAirways);
    lastID = 2; startID=lastID+1;
    save(IDlogpath,'lastID');
end

% Determine number of files to analyse
curpath = [path file(1:end-4) '/Filtered_objects/'];
folderInfo = dir(curpath); % all files in folder
nfiles = length(folderInfo);    % Number of files found
    

dispmess=msgbox({'We will now go through each image. Enter "Y" or "N" in the MATLAB window to confirm or reject objects as airways.'; ' '; 'This may take some time as there may be many objects; you will be prompted after each object to continue to the next image, or to abort the process and continue later.' ;' '; 'Please click OK to proceed.'});
drawnow; waitfor(dispmess); clear dispmess
close;

if nfiles == 2 %isnan(objvec) SHOULD BE 2 TO ACCOUNT FOR . AND ..?
    msgbox({'There are no objects to check!';'';'Run ObjectIDscript.m on your tif file to identify possible airways.'})
elseif lastID == nfiles
    error('You have checked all the objects from this tif!')
else
    disp(['There are ',num2str(nfiles-startID+1),' objects to check'])
    for objloop=startID:nfiles
        file = folderInfo(objloop).name;
        curfn=[curpath file]; %'image_Object_' num2str(obj) '.mat'];
        if exist(curfn)
            another=1;
            %disp([curfn ' exists']) DO WE NEED TO OUTPUT THIS??
            load(curfn);
            close;
            obj = file(14:end-4);
            figure('Visible','on','Name','Is this an airway?');
            imshow(zoomedimg);
            clear accept
            accept=input('Is this a suitable airway?, y or n ','s');
            % menu('Is this a suitable airway?','yes','no');
            if accept=='y'
                disp(['Object ',obj,' is an airway'])
                svpnfn=[pnAirways '/Object_', obj,'.mat'];
                save(svpnfn,'zoomedimg','lumenimg',...
                    'Blum','Llum','Nlum','Alum')
            else
                disp(['Object ',obj,' is NOT an airway'])
            end
            close all;
            clear zoomedimg lumenimg Blum Llum Nlum Alum
        else
            disp([curfn ' does not exist'])
        end
        
        % update ID log
        lastID=objloop;
        save(IDlogpath,'lastID');
        
        %Check if user wants to carry on
        disp(['There are ',num2str(nfiles-objloop),...
            ' objects left to check'])
        if mod(objloop,20)==0
            accept=histomenu('Do you want to continue?','Yes', 'No');
        end
        if accept == 2
            msgbox('OK - to continue another time from the next object, please run ConfirmAirways_script.m again')
            cd ..
            break
        end
    end
end
