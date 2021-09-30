function [tottrans,meanTheta,awdensity,flag]=TracePixelsAlongLine(...,
    pn,fn,SMAorPSR)

% [fn, pn]=uigetfile('*.mat','Select input file');

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
    hold on; plot(dat(:,2),dat(:,1),'b.'); % show start of boundary
    hold on; plot(dat(1,2),dat(1,1),'g.','MarkerSize',30);
    hold on; plot(centroid(1),centroid(2),'b*','MarkerSize',10) % plot centroid
end

% Get conversion factor
% SMAorPSR=2; % input('Is this an image of SMA or PSR?, 1=SMA, 2=PSR ');
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
            hold on; plot(xdat(itn),ydat(itn),'g.','MarkerSize',30)
        end
        itn=itn+1;
    end
    
    % Fit a line along the boundary (tangent)
    clear p x0 y0
    p=polyfit(xdat,ydat,1);
    x0=xdat(uint32(length(xdat)/2));
    y0=polyval(p,x0);
    if showfigs==1
        hold on; plot(x0,y0,'g.','MarkerSize',30)
    end
    
    % Compute and draw a line orthogonal to tangent to boundary at location
    clear pt1 pt2 v
    pt1=[xdat(1),ydat(1)];
    pt2=[xdat(end),ydat(end)];
    v=pt2-pt1; v = v / norm(v);
    if showfigs==1
        hold on; plot(linspace(xdat(1),xdat(end),100),...
            polyval(p,linspace(xdat(1),xdat(end),100)),'b-')
        hold on; line([x0+v(2), x0-v(2)],[y0-v(1), y0+v(1)]);
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
        hold on; plot(xspoke(spokenum,:),yspoke(spokenum,:),'r-')
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
                hold on; plot(xspoke(spokenum,:),...
                    yspoke(spokenum,:),'r-');
            end
        else
            % disp('First spoke is outside lumen, use it')
        end
        
        % Plot end of spoke used
        if showfigs==1
            hold on; plot(xspoke(spokenum,end),yspoke(spokenum,end),'r.',...
                'MarkerSize',10)
        end
        
        % [xi,yi]=ginput(2);
        if spokenum==1
            marker1='bo';
        else
            marker1='b*';
        end
        if showfigs==1
            hold on; plot(xspoke(spokenum,1),yspoke(spokenum,1),marker1);
            hold on; plot(xspoke(spokenum,end),yspoke(spokenum,end),'b*')
        end
        
        % Compute profiles along spokes
        [c{spokenum}]=improfile(im2bw(greycompimg,bwthresh),...
            [xspoke(spokenum,1) xspoke(spokenum,end)],...
            [yspoke(spokenum,1) yspoke(spokenum,end)]);
        
        % Compute transitions of profiles along each spoke
        %   (count changes from 0 to 1 and 1 to 0)
        transit(spokenum)=sum(diff(c{spokenum})~=0);
        % disp(['The number of transitions along spoke ',...
        %     num2str(spokenum),' is ',...
        %    num2str(transit(spokenum))])
        
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

% figure;
% for loop=1:n
%     subplot(n,1,loop)
%     pltdat=c{loop};
%     lendat=1:length(pltdat(:,1,1));
%     plot(lendat,pltdat(:,1,1),'k-');
% end

% figure; imshow(imbinarize(greycompimg,bwthresh));
% for spokenum=1:n
%     hold on;
%     plot(xspoke(spokenum,:),yspoke(spokenum,:),'r-')
%     hold on; plot(xspoke(spokenum,1),yspoke(spokenum,1),'b*');
%     hold on; plot(xspoke(spokenum,end),yspoke(spokenum,end),'b*')
% end

tottrans=sum(transit);
% disp(['Total transitions are ',num2str(tottrans)])

% disp('Angles between spokes and line from centroid to spoke starting point are ')
% disp(ThetaInDegrees')

meanTheta=mean(ThetaInDegrees);
% disp([' with average = ',num2str(meanTheta)])

end
% figure; imshow(im2bw(greycompimg,bwthresh))
% % Show lumen figure with plots
% lumfig=figure; imshow(lumenimg);
% hold on; plot(uint32(xspoke(end)),uint32(yspoke(end)),'r.','MarkerSize',10)