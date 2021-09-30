function []=Main_FromMathsServers(pn,file_tif,SMAorPSR)
close all; 
% Main function to trace SMA or PSR stains in images, using a threshold or
% a manual tracing algorithm
% ------------------------------------------------------------------------
% Custom functions called in this script:
%  * histomenu.m: function to ask for user input
%  * bezier.m: function to perform Bezier operation to alter boundaries
%  * Count_Sirius.m: function to measure area stained by PSR
%  * Count_SMA.m: function to measure area stained by alpha-SMA
%  * Threshold_LAB.m: function to threshold image using LAB mapping
%  * Threshold_greyscale.m: function to threshold image using greyscale
%  * painbrush.m: function to manually add or delete pixels from tracing
% ------------------------------------------------------------------------

% disp(['Tell me if this is an image of SMA stain or PSR stain ',...
%     'in the next dialog box.']);
% SMAorPSR=input('Is this an image of SMA or PSR?, 1=SMA, 2=PSR ');
if strcmp(SMAorPSR,'SMA') % Image is of SMA
    PixelsToMicrons=0.228; %
    % PixelsToMicrons=0.456; %
else % Image is of PSR
    PixelsToMicrons=0.228; %
end

%% File Input
results_dir = [pn '/Airways/Results'];
if isfolder(results_dir)
    rmdir(results_dir,'s')
    mkdir(results_dir)
else
    mkdir(results_dir)
end

airways_dir = [pn '/Airways/'];

folderfiles = dir([airways_dir,'Object*.mat']); % all files in folder
nfiles = length(folderfiles);    % Number of files found
if nfiles == 0
    error('There are no airways to analyse. Run ObjectIDscript.m and ConfirmAirways_script.m to identify airways.')
end

% set/load ID log for re-starting process
Airwaydata = [airways_dir,file_tif(1:length(file_tif)-4),'.mat'];
AirwayIDlogpath = [airways_dir '/AirwayIDlog.mat'];
if exist(AirwayIDlogpath)
    % load the file ID we got to last time, and update
    load(AirwayIDlogpath)
    startID=lastID+1;
    % load the data from the last run
    load(Airwaydata);
else
    lastID = 0; startID=lastID+1;
    save(AirwayIDlogpath,'lastID');
end
if lastID == nfiles
    error('You have already analysed all the available objects. If you would like to re-analyse them, simply delete the file AirwayIDlog.mat in the Airways directory.')
end

%% set up results vectors for easier writing to csv file later
% BUT DON'T NEED TO INITIALISE ANYMORE BECAUSE LOADING EXISTING DATA
% airwaywallarea = zeros(nfiles,1);
% greenarea = zeros(nfiles,1);
% innerarea = zeros(nfiles,1);
% outerarea = zeros(nfiles,1);
% lumenarea = zeros(nfiles,1); 
% epithelialarea = zeros(nfiles,1);
% R2 = zeros(nfiles,1);
% R1  = zeros(nfiles,1); 
% lumenradius = zeros(nfiles,1);
% outerperim = zeros(nfiles,1);
% innerperim = zeros(nfiles,1);
% lumenperim = zeros(nfiles,1);
% normalisedgreenarea = zeros(nfiles,1);
% normalisedepithelialarea = zeros(nfiles,1); 

%%
for file_loop=startID:nfiles
    file = folderfiles(file_loop).name;
     if startsWith(file, 'Object') 
         load([airways_dir file])
     else
         continue
     end

    close all
    
    %% now do the analysis
    scalebardisp=1; % histomenu(dispmess,'Yes','No');
    clear scalebarimg;
    scalebarimg=zoomedimg;
    convrt=PixelsToMicrons; % input('Enter conversion factor um/pixel ');
    if scalebardisp==1
        len=50; % input('Enter length of scale bar in um ');
        thick=5; % input('Enter thickness of scale bar in pixels ');
        sz=size(zoomedimg);
        x=10; % input('Enter x co-ordinate ');
        y=sz(1)-10; % input('Enter y co-ordinate ');
        x=uint32(x);y=uint32(y); % Convert x and y to real
        clear i j
        ylim=double(uint32(thick));
        xlim=double(uint32(len/convrt));
        borw=1; % histomenu('Use white or black scale bar?','White','Black');
        if borw==1
            fill=255;
        else
            fill=0;
        end
        for i=x:x+xlim
            for j=y:y+ylim
                i=uint32(i);
                j=uint32(j);
                scalebarimg(j,i,1)=fill;
                scalebarimg(j,i,2)=fill;
                scalebarimg(j,i,3)=fill;
            end
        end
    end

    %  Display image w/ scale bar
    zoomedimg=scalebarimg; clear scalebarimg;

    % Show lumen
    dispmess=msgbox('Now you will see an airway image and tracing of the lumen. When it appears use your mouse/trackpad/touchscreen to manually trace the Basement Membrane');
    drawnow; waitfor(dispmess); clear dispmess
    close;
    % Trace INNER boundary
    dat=Blum{1};
    figure;imshow(zoomedimg); hold on; 
    plot(dat(:,2),dat(:,1),'g.')
    h=drawfreehand(gca);
    clear innermask; innermask=createMask(h);
    % Find points along boundaries
    j=1; %initialize j
    i=1; %initialize i
    for i=1:length(innermask(1,:))
        for j=1:length(innermask(:,1))
            if innermask(j,i)~=0, break, end
        end
        if innermask(j,i)~=0, break, end
    end
    % Trace command
    clear innerbound; innerbound=bwtraceboundary(innermask,[j i],'N');

    close; figure('Visible','on','Name','Check to manually alter boundaries');
    imshow(zoomedimg);
    hold on
    plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);

    % Ask user to Manually alter boundaries using Bezier,
    %  first starting with inner boundary
    again=histomenu('Manually alter Epithelial Basement Membrane tracing?',...
        'yes','no');
    halfwid=50;
    while again==1
        boundary=innerbound;
        close all
        [newbound,sz,loc]=bezier(zoomedimg,boundary,halfwid);
    %     close; 
    %     figure; imshow(zoomedimg);
        figure(1), hold on
        plot(newbound(:,2),newbound(:,1),'b',...,
            'LineWidth',3,'MarkerSize',20);
        looksOK=histomenu('Keep this tracing or try again?','keep','try again');
        if looksOK==1
            boundary=newbound; clear newbound
            again=histomenu('continue to alter tracing?','yes','no');
        else
            %This uses a bezier curve to create a new boundary trace using
            %the additional point selected by the user
            changebz=histomenu('change smoothness of altered boundary tracing?','yes','no');
            if changebz==1
                disp(['Current smoothness parameter is, ',...
                    num2str(halfwid),'. Larger values increase the smoothness/width of the interpolated boundary, smaller ones result in a sharper addition to the trace.'])
                %disp(['  and length of tracing is ',num2str(sz(1)),' pixels'])
                %disp(['  and location on tracing is ',num2str(loc)])
                halfwid=input('New smoothness parameter = ');
            end
            again=1;
        end
        close all;
        innerbound=boundary;
    end

    % create new inner mask from boundary
    szinnermask=size(innermask);
    newinnermask=zeros(szinnermask(1),szinnermask(2));
    newinnermask=roipoly(newinnermask,innerbound(:,2),innerbound(:,1));
    clear innermask; innermask=newinnermask; clear newinnermask

    % dilate inner tracing to enclose airway
    outbounddilate=20*2; % outer bound dilation amount in um
    SE = strel('square',double(uint32(outbounddilate/convrt)));
    outermask=imdilate(innermask,SE);
    % Find points along boundaries
    j=1; %initialize j
    i=1; %initialize i
    for i=1:length(outermask(1,:))
        for j=1:length(outermask(:,1))
            if outermask(j,i)~=0, break, end
        end
        if outermask(j,i)~=0, break, end
    end
    % Trace boundary
    clear outerbound; outerbound=bwtraceboundary(outermask,[j i],'N');

    close;
    figure('Visible','on','Name','Check to manually alter boundaries');
    imshow(zoomedimg);
    hold on
    plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);

    plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);

    %  then alter outer boundary
    again=histomenu('Manually alter Outer Tracing of Smooth Muscle Layer?','yes','no');
    %halfwid=10;
    clear boundary newbound
    while again==1;
        boundary=outerbound;
        close all
        [newbound,sz,loc]=bezier(zoomedimg,boundary,halfwid);
    %     close;
    %     figure; imshow(zoomedimg);
        figure(1), hold on;
        plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
        hold on;
        plot(newbound(:,2),newbound(:,1),'b',...,
            'LineWidth',3,'MarkerSize',20);
        looksOK=histomenu('Keep this tracing or try again?','keep','try again');
        if looksOK==1;
            boundary=newbound; clear newbound
            again=histomenu('continue to alter tracing?','yes','no');
        else
            % Uses Bezier as above.
            changebz=histomenu('change smoothness of altered boundary tracing?','yes','no');
            if changebz==1
                disp(['Current smoothness parameter is, ',...
                    num2str(halfwid),'. Larger values increase the smoothness/width of the interpolated boundary, smaller ones result in a sharper addition to the trace.'])
                %disp(['  and length of tracing is ',num2str(sz(1)),' pixels'])
                %disp(['  and location on tracing is ',num2str(loc)])
                halfwid=input('New smoothness parameter = ');
            end
            again=1;
        end
        close all;
        outerbound=boundary;
    end
    % create new outer mask from boundary
    szoutermask=size(outermask);
    newoutermask=zeros(szoutermask(1),szoutermask(2));
    newoutermask=roipoly(newoutermask,outerbound(:,2),outerbound(:,1));
    clear outermask; outermask=newoutermask; clear newoutermask

    close;
    figure('Visible','on','Name','Final boundaries');
    imshow(zoomedimg);
    hold on;
    plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
    %hold on;
    plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);
    %hold on;
    plot(dat(:,2),dat(:,1),'g.','LineWidth',3);

    clear dilated_inner dilated_outer Bout outerboundary

    %% Run the programs to find SMA or PSR pixels and counting functions
    transl=0.6; % Set translucency of green overlay (0=transparent, 1=opaque)
    threshortrace=histomenu('Use a threshold function or manually identify areas of stain?',...
        'Threshold','Manual');
    if threshortrace==1
        % make green image
        green = cat(3, zeros(size(zoomedimg(:,:,1))), ...
            ones(size(zoomedimg(:,:,2))), zeros(size(zoomedimg(:,:,3))));
        if strcmp(SMAorPSR,'PSR') % PSR
            [level,mask]=Threshold_LAB(zoomedimg,...
                innerbound,outerbound,...
                innermask,outermask,transl,green);
            clear threshlevel; threshlevel=level;
            %disp(['Threshold levels are ',num2str(threshlevel(1)),' and ',num2str(threshlevel(2))])
        else
            [f,level(:),mask]=Threshold_greyscale(zoomedimg,...
                innerbound,outerbound,...
                innermask,outermask,transl,green);
            clear threshlevel; threshlevel=level(:);
        end
        % show user the green stuff
        clear sz; sz=size(zoomedimg);
        clear greenmask applymask displayimg
            greenmask(:,:,1)=zeros(sz(1),sz(2));
            greenmask(:,:,2)=255*mask;
            greenmask(:,:,3)=zeros(sz(1),sz(2));
            applymask(:,:,1)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,2)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,3)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask=uint8(applymask);

                close all

                figure(1), hold on
            displayimg=(zoomedimg.*applymask+uint8(greenmask));
            imshow(displayimg)
        hold on; 
        plot(innerbound(:,2),innerbound(:,1),'b','LineWidth',3);
        plot(outerbound(:,2),outerbound(:,1),'r','LineWidth',3);
        plot(dat(:,2),dat(:,1),'k','LineWidth',3);   
        alterthresh=menu('Manually alter threshold area?','yes','no');
        if alterthresh==1
            another_object=1;
            clear sz; sz=size(zoomedimg);
            clear greenmask applymask displayimg
            while another_object==1
                [mask]=paintbrush(zoomedimg,mask,transl,green);
                mask=mask.*(outermask-innermask); % only area b/w boundaries
                another_object=histomenu('Find another SMA or PSR portion?','yes','no');
            end
        end
    else
        close all
        draw_region=histomenu({'Manually outline SMA or PSR stain.';' ';'Option 1: First, use your mouse/trackpad/touchscreen to draw regions of stain';'and then edit these in detail as necessary using a "paintbrush" in the next window.';'';'Option 2: Skip drawing regions, and identify stain only with the "paintbrush"'},'Draw regions','Go straight to paintbrush');
        another_object=1;
        first_object=1;
        clear sz; sz=size(zoomedimg);
        mask=zeros(sz(1),sz(2));
        clear greenmask applymask displayimg
        while another_object==1
            if first_object==0
                draw_region=histomenu('Draw regions or go straight to editing via the "paintbrush"?','Draw regions','Paintbrush');
            end
            clear greenmask applymask displayimg
            greenmask(:,:,1)=zeros(sz(1),sz(2));
            greenmask(:,:,2)=255*mask;
            greenmask(:,:,3)=zeros(sz(1),sz(2));
            applymask(:,:,1)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,2)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,3)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask=uint8(applymask);
            displayimg=(zoomedimg.*applymask+uint8(greenmask));
            % YOU HAVE TO DRAW A REGION IN THE FIGURE THAT SHOWS UP AND
            % THEN ANOTHER GUI ALLOWS YOU TO ADD DOTS (VIA PAINTBRUSH) AND
            % THEN YOU'RE ASKED IF YOU TO DRAW ANOTHER REGION. NEED TO MAKE CLEAR TO
            % USER WHAT WILL HAPPEN ETC
            %*******            
%             figure(1), hold on
%             imshow(displayimg);
%             hold on;
%             plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
%             hold on;
%             plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);
%             clear h; h=imfreehand(gca);
%             clear tempmask; tempmask=createMask(h);
%             mask=mask+tempmask;
            % EDITS TO ALLOW FOR MULTIPLE REGIONS
            if draw_region==1
                figure(1), hold on
                imshow(displayimg);
                hold on;
                plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
                hold on;
                plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);
                clear h; h=imfreehand(gca);
                clear tempmask; tempmask=createMask(h);
                mask=mask+tempmask;
                drawanother_region=histomenu('Draw another region?','yes','no');
                while drawanother_region==1
                    clear h; h=imfreehand(gca);
                    clear tempmask; tempmask=createMask(h);
                    mask=mask+tempmask;
                    drawanother_region=histomenu('Draw another region?','yes','no');
                end
                close(1);
                dispmess=msgbox('You will now be shown the object with your selected stain region(s) highlighted. Detailed edits or further additions can be made using your mouse/touchscreen as a "paintbrush".');
                drawnow; waitfor(dispmess); clear dispmess
            end
            [mask]=paintbrush(zoomedimg,mask,transl,greenmask);
            mask=mask.*(outermask-innermask); % only area b/w boundaries
            clear greenmask applymask displayimg
            greenmask(:,:,1)=zeros(sz(1),sz(2));
            greenmask(:,:,2)=255*mask;
            greenmask(:,:,3)=zeros(sz(1),sz(2));
            applymask(:,:,1)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,2)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask(:,:,3)=ones(sz(1),sz(2)).*imcomplement(mask);
            applymask=uint8(applymask);
            displayimg=(zoomedimg.*applymask+uint8(greenmask));

            figure(1)
            imshow(displayimg)
            hold on;
            plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
            hold on;
            plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);
            another_object=histomenu('Do you want to select more SMA or PSR portions, or make further edits?','yes','no');
            if another_object==1, first_object=0;end
            close(1)
        end
    end

%     disp(['Max of mask is ',num2str(max(max(mask)))])
    for ii=1:length(mask(:,1))
        for jj=1:length(mask(1,:))
            if mask(ii,jj)>0
                mask(ii,jj)=1;
            else
                mask(ii,jj)=0;
            end
        end
    end
%     disp(['New Max of mask is ',num2str(max(max(mask)))])

    % Show and save image of tracings;
    dispmess=msgbox('A series of images will now be produced, and image data saved. Finally, you will be shown the airway object with the boundary tracings and stain areas that you have chosen.');
        drawnow; waitfor(dispmess); clear dispmess
    greenmaskimg(1:length(mask(:,1)),1:length(mask(1,:)),1)=...
        255.*(1-mask);
    greenmaskimg(1:length(mask(:,1)),1:length(mask(1,:)),2)=...
        255.*ones(length(mask(:,1)),length(mask(1,:)));
    greenmaskimg(1:length(mask(:,1)),1:length(mask(1,:)),3)=...
        255.*(1-mask);

    close all; figure; imshow(greenmaskimg);
    hold on; 
    plot(innerbound(:,2),innerbound(:,1),'b','LineWidth',3);
    plot(outerbound(:,2),outerbound(:,1),'r','LineWidth',3);
    plot(dat(:,2),dat(:,1),'k','LineWidth',3);
    clear str
    str=[results_dir,'/',file(1:length(file)-4) '_tracesonly.fig'];
    savefig(str); % save as fig
    clear str;
    str=[results_dir,'/',file(1:length(file)-4) '_tracesonly.jpg'];
    saveas(gcf,str);
    close;

    % make png
    % greenmaskpng(1:length(mask(:,1)),1:length(mask(1,:)),1)=...
    %     zeros(length(mask(:,1)),length(mask(1,:)));
    % greenmaskpng(1:length(mask(:,1)),1:length(mask(1,:)),2)=...
    %     255.*mask;
    % greenmaskpng(1:length(mask(:,1)),1:length(mask(1,:)),3)=...
    %     zeros(length(mask(:,1)),length(mask(1,:)));
    % for k=1:length(innerbound(:,1))
    %     greenmaskpng(uint32(innerbound(k,1)),uint32(innerbound(k,2)),3)=255;
    % end
    % for k=1:length(outerbound(:,1))
    %     greenmaskpng(uint32(outerbound(k,1)),uint32(outerbound(k,2)),1)=255;
    % end
    % for k=1:length(dat(:,1))
    %     greenmaskpng(uint32(dat(k,1)),uint32(dat(k,2)),3)=255;
    % end
    % figure; imshow(greenmaskpng);
    % 
    % clear str
    % str=[results_dir,'/', file(1:length(file)-4) '_tracesonly.png'];
    % greenmaskpngbw=im2bw(greenmaskpng);
    % alphamatrix=zeros(length(greenmaskpngbw(:,1)),length(greenmaskpngbw(1,:)));
    % for k=1:length(greenmaskpngbw(:,1))
    %     for j=1:length(greenmaskpngbw(1,:))
    %         if greenmaskpngbw(k,j)~=0
    %             alphamatrix(k,j)=1;
    %         end
    %     end
    % end
    % 
    % imwrite(greenmaskpng,str,'Alpha',alphamatrix)
    % close;
    %% Effectively replacing this section with OVAcompositionfunc.m 
    %
    % % Compute and display lumen Area and Perimeter
    % numberPixelsInLumen=sum(sum(lumenimg));
    % disp(['The number of pixels in the lumen is ',...
    %     num2str(numberPixelsInLumen)])
    % LumenArea=numberPixelsInLumen*(PixelsToMicrons)^2;
    % disp(['The area of the lumen is ',...
    %     num2str(LumenArea),' square microns'])
    % sz_lumenperi=size(dat);
    % LumenPerimeter=sz_lumenperi(1,1)*PixelsToMicrons;
    % disp(['The length of lumen perimeter in pixels is ',...
    %     num2str(sz_lumenperi)])
    % disp(['The length of the lumen perimeter is ',...
    %     num2str(LumenPerimeter),' square microns'])

    % now count pixels representing stains (actually only using these to get
    % final mat files etc
    if strcmp(SMAorPSR,'SMA') % Image is of SMA
        newimg = Count_SMA(zoomedimg, innermask, outermask, mask, dat, results_dir,file,PixelsToMicrons);
    else % Image is of PSR
        newimg = Count_Sirius(zoomedimg, innermask, outermask, mask, dat, results_dir,file,PixelsToMicrons);
    end

    % % Compute areas in microns^2
    % OuterArea=numberOfOuterMaskPixels*(PixelsToMicrons)^2;
    % InnerArea=numberOfInnerMaskPixels*(PixelsToMicrons)^2;

    ind = file_loop;
    filename_mat = [results_dir,'/',file(1:length(file)-4),'_final.mat'];
    filename_fig = [results_dir,'/',file(1:length(file)-4),'_tracesonly.fig'];

    [airwaywallarea(ind), greenarea(ind), innerarea(ind), outerarea(ind), lumenarea(ind), epithelialarea(ind), R2(ind), R1(ind), lumenradius(ind), outerperim(ind), innerperim(ind), lumenperim(ind), normalisedgreenarea(ind), normalisedepithelialarea(ind), ~] = OVAcompositionfunc(filename_mat,filename_fig);
    filename{ind}=file(1:length(file)-4); %FILE;

    dispmess=msgbox(['Process for airway ', num2str(file_loop), ' out of a total of ', num2str(nfiles) ' airways is complete']);
    drawnow; waitfor(dispmess); clear dispmess

    % update ID log
    lastID=file_loop;
    save(AirwayIDlogpath,'lastID');
    % save info to a mat file; construct csv table after (to allow for stop-starting of process)
    save(Airwaydata,'filename','airwaywallarea','greenarea', 'innerarea', 'outerarea', 'lumenarea', 'epithelialarea', 'R2', 'R1', 'lumenradius', 'outerperim', 'innerperim', 'lumenperim', 'normalisedgreenarea', 'normalisedepithelialarea');

    %Check if user wants to carry on
    % need to add a clause for if there are 0 objects left to check!
    carryon = histomenu(['There are ',num2str(nfiles-file_loop),' objects left to check. Do you want to continue?'],'Yes','No');

    if carryon == 2
        msgbox('OK - to continue another time from the next airway, please run AnalyseAirway_script.m again');
        cd ../
        break
    end
 
end

%% Make table and write data to csv file
% Create, populate, and save table of values

ColLabels={'filename',...
    'Airway Wall Area (um^2)',...
    'Stained Area (um^2)', ...
    'Inner Area (um^2)', ...
    'Outer Area (um^2)', ...
    'Lumen Area (um^2)', ...
    'Epithelial Area (um^2)', ...
    'Effective Outer Radius, R2 (um)',...
    'Effective Inner Radius, R1 (um)',...
    'Effective Lumen Radius, (um)', ...
    'Outer airway perimeter, (um)', ...
    'Basement Membrane perimeter (um)',...
    'Lumen perimeter (um)', ...
    'Normalised stained area (with respect to BM perimeter (um)', ...
    'Normalised epithelial area (with respect to BM perimeter (um)'};
    
clear filestr;
filestr = [airways_dir,file_tif(1:length(file_tif)-4),'.csv'];
% if exist(filestr,'file')
%     rm(filestr)
% else
    T = table(filename',airwaywallarea', greenarea', innerarea', outerarea', lumenarea', epithelialarea', R2', R1', lumenradius', outerperim', innerperim', lumenperim', normalisedgreenarea', normalisedepithelialarea');
    T.Properties.VariableNames = ColLabels;
    writetable(T,filestr,'Delimiter',',')
%end


