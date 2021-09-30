function [newimg] = Count_Sirius(img, innermask, outermask, mask, lumen,path,file,PixelsToMicrons)
% function [RegionArea,InnerPerimeter,R1,R2,Area_PSR,...
%     numberOfOuterMaskPixels,numberOfInnerMaskPixels] = ...
%     Count_Sirius(img, innermask, outermask, mask, lumen,...
%    path,file,PixelsToMicrons)

newimg=img;
sz=size(img);
PSRs=mask; % zeros(sz(1),sz(2));
for i=1:length(mask(:,1)) % length(newimg(:,1,1));
    for j=1:length(mask(1,:)) % length(newimg(1,:,1));
        if mask(i,j)==1 % fill(i,j)==1;
            % PSRs(i,j)=1;
            newimg(i,j,1)=0;
            newimg(i,j,2)=255;
            newimg(i,j,3)=0;
        else
        end
    end
end

%% Trace inner/outer boundaries

% Find points along boundaries
for i=1:length(innermask(1,:))
    % disp(['i=',num2str(i)])
    m=0; %initialize m
    for j=1:length(innermask(:,1))
        m=m+1;
        lencol(i)=m;
        if innermask(j,i)~=0, break, end  
    end
end
j=min(lencol);
i=find(lencol==j);
i=i(1); % Only use one value for starting position, i

% Trace command
innerbound=bwtraceboundary(innermask,[j i],'N');

% Plot boundaries
figure;imshow(newimg);
hold on;
plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);

% Find points along boundaries
j=1; %initialize j
i=1; %initialize i
for i=1:length(outermask(1,:))
    m=0; %initialize m
    for j=1:length(outermask(:,1))
        m=m+1;
        lencol(i)=m;
        if outermask(j,i)~=0, break, end  
    end
end
j=min(lencol);
i=find(lencol==j);
i=i(1); % Only use one value for starting position, i

% Trace boundary
outerbound=bwtraceboundary(outermask,[j i],'N');

%% Remove unwanted portions between boundaries
selimg(:,:,1)=double(img(:,:,1)).*(double(innermask==0));
selimg(:,:,1)=double(selimg(:,:,1)).*(double(outermask==1));
selimg(:,:,2)=double(img(:,:,2)).*(double(innermask==0));
selimg(:,:,2)=double(selimg(:,:,2)).*(double(outermask==1));
selimg(:,:,3)=double(img(:,:,3)).*(double(innermask==0));
selimg(:,:,3)=double(selimg(:,:,3)).*(double(outermask==1));
selimg=selimg/255;
% PSRs(:,:)=double(PSRs(:,:)).*(double(innermask==0));
% PSRs(:,:)=double(PSRs(:,:)).*(double(outermask==1));
clear newimg; newimg=img;
for i=1:length(newimg(:,1,1));
    for j=1:length(newimg(1,:,1));
        if PSRs(i,j)==1;
            newimg(i,j,1)=0;
            newimg(i,j,2)=255;
            newimg(i,j,3)=0;
        else
        end
    end
end
figure;imshow(selimg);
figure;imshow(PSRs);

%% Plot final figure
close; close; figure;imshow(newimg);
hold on;
plot(innerbound(:,2),innerbound(:,1),'g','LineWidth',3);
hold on;
plot(outerbound(:,2),outerbound(:,1),'r','LineWidth',3);
hold on;
plot(lumen(:,2),lumen(:,1),'r','LineWidth',3);

%% Set calibration factor
% PixelsToMicrons=0.228; % input('Enter number of microns in a pixel ');

%% Turn on diary and save log file
% logfilename = ... % sprintf([path file(1:length(file)-4),'_Object %d.log'], k);
%     [path,'/',file(1:length(file)-4),'.log'];
% diary(logfilename);
% diary on
% 
% %% Compute number of pixels in selected region and inner/outer perimeter
% numberOfOuterMaskPixels=sum(sum(outermask==1))
% numberOfInnerMaskPixels=sum(sum(innermask==1))
% numberPixelsInRegion=numberOfOuterMaskPixels-numberOfInnerMaskPixels;
% disp(['The number of pixels in selected region is ',...
%     num2str(numberPixelsInRegion)])
% disp(['The area of this region is ',...
%     num2str(numberPixelsInRegion*(PixelsToMicrons)^2),...
%     ' square microns'])
% sz_inner=size(innerbound);
% sz_outer=size(outerbound);
% disp(['(Manually traced) Inner perimeter is ',num2str(sz_inner(1,1)),...
%     ' pixels or ',num2str(sz_inner(1,1)*PixelsToMicrons),' microns,']);
% R1=(sz_inner(1,1)*PixelsToMicrons)/2/pi;
% disp([' so the equivalent inner radius is ',...
%     num2str(R1),' microns']);
% disp(['(Manually traced) Outer perimeter is ',num2str(sz_outer(1,1)),...
%     ' pixels or ',num2str(sz_outer(1,1)*PixelsToMicrons),' microns,']);
% R2=(sz_outer(1,1)*PixelsToMicrons)/2/pi;
% disp([' so the equivalent outer radius is ',...
%     num2str(R2),' microns']);
% PSR_pixels=sum(sum(PSRs));
% disp(['Total number of pixels representing PSRs is ',num2str(PSR_pixels)])
% Area_PSR=PSR_pixels*(PixelsToMicrons)^2; %microns^2
% disp(['So the area of airway stained by picrosirius red is ',...
%     num2str(Area_PSR),' microns^2']);
% RegionArea=numberPixelsInRegion*(PixelsToMicrons)^2;
% InnerPerimeter=sz_inner(1,1)*PixelsToMicrons;
%     
% %% Turn off diary
% diary off

%% Save figures and data to file
sv=menu('save figure and data?','yes','no');
if sv==1
    % Save final image
    clear str
    str=[path,'/',file(1:length(file)-4),'_final.mat'];
    save(str,'newimg')
    
    close all;
    % Raw Image
    filename = [path,'/',file(1:length(file)-4),'.tif'];
    imwrite(img,filename,'tif');
      
    % Traces of Boundaries
    figure;imshow(newimg);
    hold on;
    plot(innerbound(:,2),innerbound(:,1),'b','LineWidth',3);
    hold on;
    plot(outerbound(:,2),outerbound(:,1),'r','LineWidth',3);
    hold on;
    plot(lumen(:,2),lumen(:,1),'k','LineWidth',3);
    clear str;
    str=[path,'/',file(1:length(file)-4),'_traces.fig'];
    savefig(str);
    clear str;
    str=[path,'/',file(1:length(file)-4) '_traces.jpg'];
    saveas(gcf,str);
    close;
end

end