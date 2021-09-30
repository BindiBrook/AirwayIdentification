% code to test global vs local scaling and maybe also Mike's vs Sam's code

%% Choose the minima and maxima for area and perimeter

magnification=40;
PixelsToMicrons=0.228;
ratio_lowerlimit=40;
    
% Set minimum and maximum allowable object (airway) sizes
MinEffDia=25; % 50; % % set minimum effective diameter
MaxEffDia=500; % 275; % 250; % set maximum effective diameter

% Display to user
disp(['Range of effective diameters is [',num2str(MinEffDia),...
    ',',num2str(MaxEffDia),'] microns'])

% Compute Min,Max Area in Microns for given effective diameters
MinAreaMicrons=pi*(MinEffDia/2)^2;
MaxAreaMicrons=pi*(MaxEffDia/2)^2;

disp(['Range of Areas is [',num2str(MinAreaMicrons),...
    ',',num2str(MaxAreaMicrons),'] square microns'])

% Convert to min,max areas of airways to include as objects
minarea=MinAreaMicrons/((PixelsToMicrons)^2);
maxarea=MaxAreaMicrons/((PixelsToMicrons)^2);

%% obtain global scaling 
I=imread('D34_5.2_15-2018-03-09_aSMA_tile_97_rgb.tif');
grey_unsc = rgb2gray(I);
global_lims = stretchlim(grey_unsc)

imgsegsz = size(I);
w = imgsegsz(2);
h = imgsegsz(1);

%% obtain scaling based on shrunk image (called it local but I don't really mean that)
fun = @(block_struct) imresize(block_struct.data,0.01);
I2 = blockproc('D34_5.2_15-2018-03-09_aSMA_tile_97_rgb.tif',[5000 5000],fun);
grey = rgb2gray(I2);
local_lims = stretchlim(grey)


%% apply local and global scaling?
gray_global = imadjust(grey_unsc, global_lims, []);
gray_local = imadjust(grey_unsc, local_lims, []);

figure
imshow(gray_global)
figure
imshow(gray_local)

%% binarize the images and compare

BW_global = imbinarize(gray_global, 'adaptive', 'Sensitivity', 0.88);
BW_local = imbinarize(gray_local, 'adaptive', 'Sensitivity', 0.88);

figure
imshow(BW_global)
figure
imshow(BW_local)

%% So these next two imerode calls can be modified to improve/worsen results
% You can change the size of the structuring element from 2 up to
% anything,I'd stay in the 2-10 range though....
% You can also add another parameter to the imerode() call, an integer comma
% separated after the 'se' that tells the imerode function how many times to
% run.  I'm running mine twice here, and with two different calls, because I
% was visualising the effects.
% for kk=2:10;
BW = BW_global;

kk = 7;
se = strel('disk',kk);
erodedimg = imerode(BW, se); % clear BW; % close breaks in the epithelium
% eimgfname = [output_dir file(1:10) '/' file(1:length(file)-4),...
%     '_erodedimg1_se', num2str(kk), '.tif'];
% imwrite(erodedimg, eimgfname);

erodedimg2 = erodedimg;
% erodedimg2 = imerode(erodedimg, se); % close breaks in the epithelium
% clear erodedimg;
% eimg2fname = [output_dir file(1:10) '/' file(1:length(file)-4),...
%     '_erodedimg2_se', num2str(kk), '.tif'];
% imwrite(erodedimg2, eimg2fname);

% end


%This is one of those times I had to invert the image
filledborders = imfill(imcomplement(erodedimg2), [1,1;h,w;h,1;1,w]);
filledborders = imcomplement(filledborders);
% fnametemp = [output_dir file(1:10),'/','block_' num2str(location(1)),...
% '__', num2str(location(2)),'_g_fill.tif'];
% imwrite(filledborders, fnametemp);

openimg = bwareaopen(filledborders, 1000);

% This morph runs three times, as dilate's single run effects are pretty small
dilatedimg = bwmorph(openimg, 'dilate', 3);

% fname7 = [output_dir file(1:10),'/','block_' num2str(location(1)), '__', num2str(location(2)),'_i_dilate.tif'];
% imwrite(dilatedimg1, fname7);

% Clears out the stray pixels in the middle of the white objects.
processedimg = imfill(dilatedimg,'holes');

figure
imshow(processedimg)

%% Code where objects outside range [minarea maxarea] are omitted
[label, NumConnectObjs] = bwlabel(processedimg); % label the processed image

props = regionprops(logical(processedimg), ...
                    processedimg, 'Area'); % get region properties

areas = [props.Area]; % get Area

% Select objects within the range [minarea maxarea]
allowableAreaIndices = (areas>minarea) & (areas<maxarea);

% Give the indices for objects matching area criteria
keeperIndices = find(allowableAreaIndices);

% Below is the new image in area range
keeperBlobsImage = ismember(label, keeperIndices);

% Finds the boundaries of the found objects, white
[B, labelledobjects] = bwboundaries(keeperBlobsImage, 'noholes');
% logic_labelobjects=logical(keeperBlobsImage);
    % MAY WANT TO CHANGE THE ABOVE IF CURRENT RESULTS DONT MATCH
    % ANDREW'S RESULTS

% Finds information about all of the bounded objects in the image
stats = regionprops(labelledobjects,...
                    'Centroid', 'Area', 'Perimeter', 'Extrema', ...
                    'MajorAxisLength', 'MinorAxisLength', 'Orientation');

% Get number of objects
NumberofObjects = length(B);

disp(['We have identified ', num2str(NumberofObjects), ...
      ' objects that might be airways in this tile of the image ']);
disp(['   in effective diameter range: [', num2str(MinEffDia), ...
      ',', num2str(MaxEffDia), '] microns, out of the ', ...
      num2str(NumConnectObjs), ' connected objects (total blobs)']);
% disp('Press any button to continue');
% pause
% disp('registered button press, thanks!')
disp(['Now using automated process to search the ', ...
    num2str(NumberofObjects), ...
    ' identified objects for those most likely ']);
disp(' to be airways, based on circularity parameter');
% disp('this could take quite some time if there are many objects')

% Creates an invisible figure then plots x's over the spots found
% Then saves the figure to a file to compare with the original image
fig = figure('visible','off');
fig.Units = 'normalized';
fig.Position = [0 0 1 1];
redamt = 1/8;
% scaledimage = imresize(image,redamt);
% close all; 
% curfig = figure; 
% imshow(scaledimage);
% 
% hold on

% Initialise the InterestingDataArray
InterestingDataArray=[];
