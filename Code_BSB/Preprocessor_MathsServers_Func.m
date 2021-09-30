function returnVal = Preprocessor_MathsServers_Func(img_full_path,output_dir,SMAorPSR,tile_size, overlap)
% Preprocessor
%   This is the preprocessor used to identify potential airways.
%   It locates potential airways and saves results to a file
%   Then, the user runs the "Main.m" program to analyse images.

% objectCount is used to keep a running total of objects found over all
% tiles - this is useful for unique naming of files across tiles.
global objectCount
objectCount = 0;
global tileCount
tileCount = 0;
disp('Begin Preprocessor_MathsServers_Func')
% Make sub-directory called "objects"
mkdir([output_dir '/objects/']);
[~, file, ~] = fileparts(img_full_path);

% Turn on diary to record output and parameters
diary([output_dir '/' file, '.log']);
diary on

% Display time at start of program
starttime=clock;
disp('Starting program at')
disp(['Year = ',num2str(starttime(1)),...
    ', Month = ',num2str(starttime(2)),...
    ', Day = ',num2str(starttime(3)),...
    ', Hour = ',num2str(starttime(4)),...
    ', Minute = ',num2str(starttime(5)),...
    ', Second = ',num2str(starttime(6))])

tic
% img_full_path=[thispath,'/', file];

disp(['Image file path is ',img_full_path]);


%% Get image magnification and conversion factors
if strcmp(SMAorPSR,'SMA') % Image is of SMA
    magnification=40;
    PixelsToMicrons=0.228;
    ratio_lowerlimit=40;
elseif strcmp(SMAorPSR,'PSR') % Image is of PSR
    magnification=40;
    PixelsToMicrons=0.228;
    ratio_lowerlimit=40; % 60; % 80/2;
else
    disp("ERROR!  Specify if image is SMA or PSR using the input 'SMAorPSR'");
    return
end
disp(['Andrews ratio lower limit is ',num2str(ratio_lowerlimit)])


%% Choose the minima and maxima for area and perimeter
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

disp(['Range of areas in pixels is [',num2str(minarea),...
    ',',num2str(maxarea),'] pixels'])

% Set overlap amount in microns, unless it's already been supplied to the
% function
if isequal(overlap, [-1, -1])
%    overlap_scalar = double(ceil(500/PixelsToMicrons)); % SET OVERLAP TO 500 microns
    overlap_scalar = double(ceil(100/PixelsToMicrons)); % SET OVERLAP TO 100 microns
    overlap = [overlap_scalar overlap_scalar];
end
% NEED CODE IN SELECT AIRWAYS PROGRAM TO ELIMINATE DUPLICATES!


%% Declare the `preproc` function handle, which points to preproc_inner.
preproc = @(block_struct) preproc_inner(block_struct.data, block_struct.location, block_struct.imageSize, block_struct.blockSize, block_struct.border, output_dir, file, minarea, maxarea, MinEffDia, MaxEffDia, magnification, ratio_lowerlimit);
imginfo = imfinfo(img_full_path);
imginfo = imginfo(1);

if imginfo.Width < tile_size(2)
    disp(['tile_size ',num2str(tile_size(2)),' is larger than info width ',num2str(imginfo.Width)])
    disp('Aborting, because otherwise we get many duplicates.')
    returnVal = 1;
    return
end
if imginfo.Height < tile_size(1)
    disp(['tile_size ',num2str(tile_size(1)),' is larger than info height ',num2str(imginfo.Height)])
    disp('Aborting, because otherwise we get many duplicates.')
    returnVal = 1;
    return
end

% Apply preproc() to the image using blockproc.
blockproc(img_full_path, tile_size, preproc, "BorderSize", overlap);


%% Combine each of the InfoFile_*.mat files (which each correspond to a tile
% from the blockproc) into one InfoFile.mat.
file_list = dir([output_dir '/InfoFile_*.mat']);
filenames = {file_list(:).name};
InterestingDataArray = [];
for i = 1:numel(filenames)
    newData1 = load([output_dir '/' filenames{i}]);
    InterestingDataArray = cat(1, InterestingDataArray, newData1.InterestingDataArray);
end

save([output_dir '/InfoFile.mat'],...
    'file', 'InterestingDataArray');


% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end


disp(['Image was split into ',num2str(tileCount),' tiles'])
disp(['Overall, ',num2str(objectCount),' (possibly duplicate) objects were detected.'])
%% Display time program ended
endtime=clock;
disp('End is at ')
disp(['Year = ',num2str(endtime(1)),...
    ', Month = ',num2str(endtime(2)),...
    ', Day = ',num2str(endtime(3)),...
    ', Hour = ',num2str(endtime(4)),...
    ', Minute = ',num2str(endtime(5)),...
    ' Second = ',num2str(endtime(6))])

disp('DONE. Now check and sort duplicates with CheckForDuplicates.m, ')
disp(' then filter potential airways with FilterObjectFiles.m, ')
disp([' and finally, SelectAirways_FromMathsServers.m, ',...
            'then Main_FromMathsServers.m'])
diary off
returnVal = 0;
return

