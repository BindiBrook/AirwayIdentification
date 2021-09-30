function [airwaywallarea, greenarea, innerarea, outerarea, lumenarea, epithelialarea, R2, R1, lumenradius, outerperim, innerperim, lumenperim, normalisedgreenarea, normalisedepithelialarea, fractiongreenpixels] = OVAcompositionfunc(matfile,figfile)
%-----------------
% Inputs:
% matfile               - _final.mat file (points to relevant directory)
% figfile               - _tracesonly.fig file (points to relevant directory)

% Outputs
% lumen_radius          - Effective lumen radius
% fractiongreenpixels   - area fraction of ASM or PSR


load(matfile) %load image file

% scaling to convert pixels to lengths (microns)
pixelstomicrons=0.228;

%% pick out the green pixels

% indices of green pixels (SMA or PSR)
notred=newimg(:,:,1)==0;
notblue=newimg(:,:,3)==0;
[iind,jind] = find(notred.*notblue);
mask = (notred.*notblue)>0;

% coordinates of green pixels
xcoord = pixelstomicrons*iind;
ycoord = pixelstomicrons*jind;

%% load up fig file to get boundary; use to construct annuli to count pixels
fig = openfig(figfile);

ax=fig.Children;
ln=ax.Children;

xin=ln(3).XData;  %get inner boundary coordinates from figure 
yin=ln(3).YData;
xlum=ln(1).XData;  %get lumen boundary coordinates from figure 
ylum=ln(1).YData;
gcf; close

% set up annulus and inner boundary mask
outbound=linspace(0,40,21);
step = outbound(2)-outbound(1); if ~(round(step)-step)==0, error('STEP SHOULD BE INTEGER'),end
convrt=pixelstomicrons;

innermask=roipoly(newimg,xin,yin); % sets 1s for pixels inside the basement membrane trace, 0s outside
SE = strel('square',double(uint32(outbound(end)/convrt)));
outermask=imdilate(innermask,SE);

% get lumen mask
lumenmask=roipoly(newimg,xlum,ylum); % sets 1s for pixels inside the lumen trace, 0s outside

%% loop over this dilation to get distribution as function of r
for i=1:length(outbound)-1
    outbounddilate=outbound(i); % outer bound dilation amount in um
    SE1 = strel('square',double(uint32(outbounddilate/convrt)));
    SE2 = strel('square',double(uint32((outbounddilate+step)/convrt)));

    if i==1, innermasknew=innermask; else, innermasknew = imdilate(innermask,SE1); end
    outermasknew = imdilate(innermask,SE2);
    annulus = outermasknew-innermasknew;

    numberOfOuterMaskPixels(i)=sum(sum(outermasknew==1));
    numberOfInnerMaskPixels(i)=sum(sum(innermasknew==1));
    numberPixelsInRegion(i)=numberOfOuterMaskPixels(i)-numberOfInnerMaskPixels(i);

    numberofgreenpixels(i)=sum(sum(mask.*annulus));

end
fractiongreenpixels=numberofgreenpixels./numberPixelsInRegion;

%% compute airway data for excel
%area
lumenpixels = sum(sum(lumenmask==1));
totaloutermaskpixels = sum(sum(outermask==1));
totalinnermaskpixels = sum(sum(innermask==1));
totalairwaypixels = totaloutermaskpixels-totalinnermaskpixels;
airwaywallarea=totalairwaypixels*convrt^2; %wall area between BM and dilated BM
greenarea = sum(numberofgreenpixels)*convrt^2;  %area of stained pixels (SMA or PSR)
innerarea = totalinnermaskpixels*convrt^2; %area inside BM
outerarea = totaloutermaskpixels*convrt^2;  % area inside dilated BM
lumenarea = lumenpixels*convrt^2;   % area of lumen
epithelialarea=innerarea-lumenarea; % area inside BM - lumen area to give epithelial cell area

%radii
R2 = sqrt(totaloutermaskpixels*convrt^2/pi); % radius of dilated BM based on area
R1 = sqrt(totalinnermaskpixels*convrt^2/pi);  % radius of BM based on area
lumenradius = sqrt(lumenpixels*convrt^2/pi); % radius of lumen based on area

%via boundary "length"
outerbound = boundarylength(outermask);
innerbound = boundarylength(innermask);
lumenbound = boundarylength(lumenmask);

%perimeters
outerperim = boundaryperimeter(outerbound,convrt); 
innerperim = boundaryperimeter(innerbound,convrt);
lumenperim = boundaryperimeter(lumenbound,convrt);

normalisedgreenarea=greenarea/innerperim; %stained area normalised to BM perimeter
normalisedepithelialarea=epithelialarea/innerperim; %epithelial area normalised to BM perimeter
