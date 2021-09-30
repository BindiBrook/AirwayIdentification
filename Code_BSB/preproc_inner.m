function preproc_inner(image, block_location, imageSize, blockSize, border, output_dir, file, minarea, maxarea, MinEffDia, MaxEffDia, magnification, ratio_lowerlimit)
% objectCount is a global variable used to track the number of objects over
% all tiles. This logic will likely break if used in parallel.
global objectCount
global tileCount
tileCount = tileCount + 1;

%% Calculate the cropping required to remove any pixels beyond the true
% image boundary.
% Also calculate the offset required to apply to object coordinates due to
% this transformation to get correct global coordinates.
remove_top = 0;
remove_bottom = 0;
remove_left = 0;
remove_right = 0;

if block_location(1) == 1
    remove_top = border(1);
    offset(1) = 0;
else
    offset(1) = border(1);
end

if block_location(2) == 1
    remove_left = border(2);
    offset(2) = 0;
else
    offset(2) = border(2);
end

if block_location(1) + blockSize(1) + border(1) > imageSize(1)
    remove_bottom = block_location(1) + blockSize(1) + border(1) - imageSize(1) - 1;
end

if block_location(2) + blockSize(2) + border(2) > imageSize(2)
    remove_right = block_location(2) + blockSize(2) + border(2) - imageSize(2) - 1;
end

%% Apply cropping to edges of image
image = image(1+remove_top:end-remove_bottom, 1+remove_left:end-remove_right, 1:3);


%% Get image info (width & height)
imgsegsz = size(image);
w = imgsegsz(2);
h = imgsegsz(1);

tic

%% Save rgb tifs
% fnamergb = [output_dir '/' file '_tile_' num2str(tileCount), '_rgb.tif'];
% imwrite(image, fnamergb);

%% Convert to greyscale and binarize
grey = rgb2gray(image);
%fname1 = [output_dir '/' file '_tile_' num2str(tileCount), '_grey1.tif'];
%imwrite(grey, fname1);
toc
% disp('Image converted to greyscale')
% clear image
% disp('Image converted to grayscale and cleared, Press any button to continue');
% pause
% disp('registered button press, thanks!')

tic
gray2 = imadjust(grey, stretchlim(grey), []);
%stretchlim(grey)
%gray2 = imadjust(grey, [0.6471 0.851], []);  %testing global scaling
%gray2 = imadjust(grey, [0.4118 0.9057], []);  %testing global scaling
%gray2 = imadjust(grey, [0.43 0.87], []);  %testing global scaling

%fname2 = [output_dir '/' file '_tile_' num2str(tileCount), '_grey2.tif'];
%imwrite(gray2, fname2);
toc
% disp('finished imadjust')

tic
BW = imbinarize(gray2, 'adaptive', 'Sensitivity', 0.88);
% fname3 = [output_dir file(1:10), '/','block_' num2str(location(1)), '__', num2str(location(2)),'_d_BW.tif'];
% imwrite(BW, fname3);

%% So these next two imerode calls can be modified to improve/worsen results
% You can change the size of the structuring element from 2 up to
% anything,I'd stay in the 2-10 range though....
% You can also add another parameter to the imerode() call, an integer comma
% separated after the 'se' that tells the imerode function how many times to
% run.  I'm running mine twice here, and with two different calls, because I
% was visualising the effects.
% for kk=2:10;
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
toc

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
procimgfname = [output_dir '/' ...
                file '_tile_' num2str(tileCount), ...
                '_procimg.tif'];
%imwrite(processedimg, procimgfname);

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
toc

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
scaledimage = imresize(image,redamt);
close all; 
curfig = figure; 
imshow(scaledimage);

hold on

% Initialise the InterestingDataArray
InterestingDataArray=[];

%% Loop over all objects
for k = 1:NumberofObjects

    % countdown=length(B)-k;
    % disp(['Analysed object ',num2str(k),...
    %     ' now there are ',num2str(countdown),' left to go']);

    % obtain (X,Y) boundary coordinates corresponding to spot 'k'
    boundary = B{k};

    % Obtain labelled airways from boundary trace
    Lumen = (labelledobjects==k);

    % compute a simple estimate of the object's perimeter
    % Remnant... RegionProps has perimeter estimates..See P2
    % delta_sq = diff(boundary).^2;
    % perimeter = sum(sqrt(sum(delta_sq,2)));
    p2 = stats(k).Perimeter;

    % obtain the area calculation corresponding to label 'k'
    area = stats(k).Area;

    % compute metrics for identifying airways
    % This is a new metric based on major/minor axes ratios
    MajorAxisLength = stats(k).MajorAxisLength; 
    MinorAxisLength = stats(k).MinorAxisLength;
    axismetric = abs(1-MajorAxisLength/MinorAxisLength);
    % This is Lesicko's circularity metric
    % metric = 4*pi*area/p2^2;
    % This metric is Andrew's metric
    areaperimratio = area/p2; % area/perimeter;
    % Perimeter is 2*pi*r, so diamater is Perimer/pi
    diameter = floor(p2/pi);

    %% Select objects as airways
    if areaperimratio > ratio_lowerlimit % && axismetric < threshold % || smetric > apratio
        % Dilate the lumens
        dilatedimg = bwmorph(Lumen, 'dilate');
        
        % Outer Mask
        M = 300*(magnification/40);
        se2 = ones(2*M + 1, 2*M + 1);
        dilated_outer = imdilate(dilatedimg,se2);

        disp(['The axismetric for object ', num2str(k+objectCount), ...
              ' is ', num2str(axismetric)])
        disp(['   and the areaperimratio for object ', num2str(k+objectCount), ...
              ' is ', num2str(areaperimratio)])
        
        centroid = stats(k).Centroid;
        c1 = floor(centroid(1));
        c2 = floor(centroid(2));
        
        % Calculate initial boundary
        [Bout, Lout, Nout, Aout] = bwboundaries(dilated_outer,'noholes');
        outerboundary = Bout{1};
        % Zoom in on kth object
        zoomedimg = image(min(outerboundary(:,1)):max(outerboundary(:,1)), ...
                          min(outerboundary(:,2)):max(outerboundary(:,2)), :);
        lumenimg = Lumen(min(outerboundary(:,1)):max(outerboundary(:,1)), ...
                         min(outerboundary(:,2)):max(outerboundary(:,2)), :);

        % Trace the lumen
        [Blum, Llum, Nlum, Alum] = bwboundaries(lumenimg, 'noholes');

        %% Add data about this object to InterestingDataArray
        addrow=[0, 0, k+objectCount, c2 + block_location(1) - 1 - offset(1), c1 + block_location(2) - 1 - offset(2), diameter];
        InterestingDataArray = [InterestingDataArray; addrow];
        
        %% Save data about this object into its own .mat file
        svpnfn=[output_dir '/objects/image_Object_' num2str(k+objectCount) '.mat'];
        
        save(svpnfn, 'zoomedimg', 'lumenimg', ...
             'Blum', 'Llum', 'Nlum', 'Alum', 'block_location')
         
        % Plot locations of objects (potential airways)
        % Add text to image
        scaledlabels = imresize(labelledobjects==k, redamt);
        txt=['Object ', num2str(k+objectCount)];
        [row,col] = find(scaledlabels==1);
        text(col(1), row(1), txt);
        hold on;
    end   
end

%% Increment global objectCount to give unique numbering of objects (so long
% as we run in serial)
objectCount = objectCount + NumberofObjects;

%% Save InterestingDataArray
save([output_dir '/InfoFile_' num2str(block_location(1)) '_' num2str(block_location(2)) '.mat'], ...
     'file', 'InterestingDataArray');
%% Save labelled objects figure

labelimgfname = [output_dir '/' file '_tile_' num2str(tileCount), '_LabelledObjects.fig'];
savefig(curfig,labelimgfname);
end