function []=FilterObjectFiles(pn,SMAorPSR)

warning('off','all'); % turn off warnings;

%% Set ranges
%ranges from Mike's method's write-up (but density lower limit changed from 80 to 75
if SMAorPSR==1
    trans_low=0;
    trans_up=501;
    theta_low=0;
    theta_up=111;
    dens_low_lib=75;  % this was 80
    dens_up_lib=125;
    dens_low_cons=75; % this was 80
    dens_up_cons=125;  
else
% ranges from Mike's original coode
%if SMAorPSR==2  %% PSR 
    trans_low=0;
    trans_up=90;
    theta_low=0;
    theta_up=90;
    dens_low_lib=83; %95;
    dens_up_lib=141; %125;  Making these the same as cons values to do away with Folder 1 and Folder 2
    dens_low_cons=83;
    dens_up_cons=141;
end

objects_dir = [pn '/objects/'];

% make subdirectories
mkdir([pn '/Filtered_objects']);
%mkdir([pn '/Folder1']);
mkdir([pn '/Folder2']);

disp('pn value')
disp(pn)
disp('SMAorPSR value')
disp(SMAorPSR)

% Get all files in this directory
clear folderfiles nfiles
folderfiles = dir(objects_dir); % all files in folder

nfiles = length(folderfiles);    % Number of files found

for ii=1:nfiles
    clear currentfilename tottrans meanTheta pxdensity
    currentfilename = folderfiles(ii).name;
    if startsWith(currentfilename, 'image') && currentfilename(end-4)~='e'
        [tottrans,meanTheta,pxdensity]= TracePixelsAlongLine(objects_dir,currentfilename,SMAorPSR); 
            
            
        if isnan(tottrans)~=1 && isnan(meanTheta)~=1 &&  isnan(pxdensity)~=1
                 
            if tottrans>=trans_low && tottrans<=trans_up &&  meanTheta>=theta_low && meanTheta<=theta_up && pxdensity>=dens_low_cons && pxdensity<=dens_up_cons
                         
                load([objects_dir currentfilename])
                if pxdensity>=dens_low_lib && pxdensity<=dens_up_lib
                    % place the file in folder 1
                    copyfile([objects_dir currentfilename],[pn '/Filtered_objects/'])
                       
                    disp([objects_dir currentfilename ' into the "Filtered_objects" folder'])
                        
                else
                    % place the file in folder 2
                    copyfile([objects_dir currentfilename],[pn '/Folder2/'])
                        
                    disp([objects_dir currentfilename ' into Folder2'])
                        
                end
            end
        end
    end
end

end


function [tottrans,meanTheta,awdensity,flag]=TracePixelsAlongLine(pn,fn,SMAorPSR)
    
load([pn fn]); % load file

showfigs=0; % turn on/off figure display

flag=0; % initialise flag to zero

n=15; % number of spokes pointing outward from lumen
t=linspace(0,1,n+1); % evenly spaced positions on [0,1] interval
t=t(1:(end-1)); % eliminate start and end points

nboundpxls=100; % set number of pixels in boundary region for computing spoke
% Set size of spoke
olen=100; % length fo spoke
ndat=100; % number of pixels making up spoke (for plotting only)

% Convert image complement to grayscale
greycompimg=imcomplement(rgb2gray(zoomedimg)); % edge(rgb2gray(zoomedimg),'Sobel'); %

% Find centroid of lumen
stats = regionprops(lumenimg);
centroid = stats.Centroid;

% Get lumen boundary
dat=Blum{1}; % with lumen boundary tracing

% Show image
if showfigs==1
    figure; imshow(greycompimg);
    plot(dat(:,2),dat(:,1),'b.'); % show start of boundary
    plot(dat(1,2),dat(1,1),'g.','MarkerSize',30);
    plot(centroid(1),centroid(2),'b*','MarkerSize',10) % plot centroid
end

% Get conversion factor
if SMAorPSR==1 % Image is of SMA
    % PixelsToMicrons=0.456; %
    PixelsToMicrons=0.228; %
else % Image is of PSR
    PixelsToMicrons=0.228; %
end
convrt=PixelsToMicrons; % conversion factor um/pixel

% Set black/white threshold value
bwthresh=50/256;

xdat=ones(nboundpxls,1); ydat=xdat; % initialize vectors
clear theta; theta=ones(1,n); % initialize vectors
for spokenum=1:n % iterate through spokes
    if spokenum==1
        spokepos(spokenum)=1;
    else
        spokepos(spokenum)=uint32(t(spokenum)*length(dat));
    end
    
    % Select some pixels along boundary
    itn=1; clear xdat ydat
    for k=spokepos(spokenum):spokepos(spokenum)+nboundpxls
        if k>length(dat(:,1))
            flag=1;
            tottrans=NaN;
            meanTheta=NaN;
            awdensity=NaN;
            return
        end
        xdat(itn)=dat(k,2); ydat(itn)=dat(k,1);
        if showfigs==1
           plot(xdat(itn),ydat(itn),'g.','MarkerSize',30)
        end
        itn=itn+1;
    end
    
    % Fit a line along the boundary (tangent)
    clear p x0 y0
    p=polyfit(xdat,ydat,1);
    x0=xdat(uint32(length(xdat)/2));
    y0=polyval(p,x0);
    if showfigs==1
        plot(x0,y0,'g.','MarkerSize',30)
    end
    
    % Compute and draw a line orthogonal to tangent to boundary at location
    clear pt1 pt2 v
    pt1=[xdat(1),ydat(1)];
    pt2=[xdat(end),ydat(end)];
    v=pt2-pt1; v = v / norm(v);
    if showfigs==1
         plot(linspace(xdat(1),xdat(end),100),...
            polyval(p,linspace(xdat(1),xdat(end),100)),'b-')
        line([x0+v(2), x0-v(2)],[y0-v(1), y0+v(1)]);
    end
    
    % Compute angle between orthogonal line and x-axis
    theta(spokenum)=(atan(v(1)/(v(2))));
    
    % Create spoke y coordinates
    clear pspoke
    pspoke=polyfit(linspace(x0-v(2),x0+v(2),ndat),...
        linspace(y0+v(1),y0-v(1),ndat),1);
    
    % First spoke x coords
    % compute spoke coordinates
    xspoke(spokenum,:)=linspace(x0,x0+olen*cos(theta(spokenum)),ndat);
    yspoke(spokenum,:)=polyval(pspoke,xspoke(spokenum,:));
    if showfigs==1
        plot(xspoke(spokenum,:),yspoke(spokenum,:),'r-')
        %hold on;
    end
    
    % Check if end of spoke is inside or outside lumen
    if uint32(yspoke(spokenum,end))>=length(lumenimg(:,1)) || ...
        uint32(xspoke(spokenum,end))>=length(lumenimg(1,:)) || ...
        uint32(yspoke(spokenum,end))<1 || uint32(xspoke(spokenum,end))<1
        flag=1;
        tottrans=NaN;
        meanTheta=NaN;
        awdensity=NaN;
        return
    else
        if lumenimg(uint32(yspoke(spokenum,end)),...
                uint32(xspoke(spokenum,end)))==1
            % spoke is in lumen, so don't use first spoke
            % disp('First spoke is in lumen, so use second spoke')
            % Second spoke x coords
            % compute spoke coordinates
            xspoke(spokenum,:)=linspace(x0,x0-olen*cos(theta(spokenum)),ndat);
            yspoke(spokenum,:)=polyval(pspoke,xspoke(spokenum,:));
            if showfigs==1
                 plot(xspoke(spokenum,:),...
                    yspoke(spokenum,:),'r-');
            end
        else
            % disp('First spoke is outside lumen, use it')
        end
        
        % Plot end of spoke used
        if showfigs==1
            plot(xspoke(spokenum,end),yspoke(spokenum,end),'r.',...
                'MarkerSize',10)
        end
        
        if spokenum==1
            marker1='bo';
        else
            marker1='b*';
        end
        if showfigs==1
            plot(xspoke(spokenum,1),yspoke(spokenum,1),marker1);
            plot(xspoke(spokenum,end),yspoke(spokenum,end),'b*')
        end
        
        % Compute profiles along spokes
        [c{spokenum}]=improfile(im2bw(greycompimg,bwthresh),...
            [xspoke(spokenum,1) xspoke(spokenum,end)],...
            [yspoke(spokenum,1) yspoke(spokenum,end)]);
        
        % Compute transitions of profiles along each spoke
        %   (count changes from 0 to 1 and 1 to 0)
        transit(spokenum)=sum(diff(c{spokenum})~=0);
        
        % Compute angle between spoke and line from centroid to start of spoke
        u1(spokenum,:)=[(centroid(1)-xspoke(spokenum,1)) ...
            (centroid(2)-yspoke(spokenum,1)) 0];
        u2(spokenum,:)=[(xspoke(spokenum,1)-xspoke(spokenum,end)) ...
            (yspoke(spokenum,1)-yspoke(spokenum,end)) 0];
        ThetaInDegrees(spokenum) = atan2d(norm(cross(u1(spokenum,:),u2(spokenum,:))),...
            dot(u1(spokenum,:),u2(spokenum,:)));
    end
end

% dilate lumen tracing to enclose airway wall
outbounddilate=20*2; % outer bound dilation amount in um
SE = strel('square',double(uint32(outbounddilate/convrt)));

% get mask for region including airway wall and lumen
outermask=imdilate(lumenimg,SE); 

% Get mask for only airway wall
mask=outermask-lumenimg;

% get only airway wall using masks
airwaywall=uint8(mask).*greycompimg; % uint8(imbinarize(greycompimg,bwthresh)); 
    % this is the region including airway wall

% Compute pixel density in airway wall
awdensity=sum(sum(double(airwaywall)))/sum(sum(double(mask)));

tottrans=sum(transit);

meanTheta=mean(ThetaInDegrees);

end
