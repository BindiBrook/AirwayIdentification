function []=CheckForDuplicates(pn,fn)
% clearvars; close all; clc

%% File Input
load([pn '/' fn]);

disp(['Running airway set: ' pn])

% Go through file and display possible duplicates
itn=0;
nfiles=size(InterestingDataArray,1);
orig = [];
for ii=1:nfiles
    currentobject=InterestingDataArray(ii,3); % get current object for comparison
    for jj=1:nfiles % iterate through files in folder
        evalobj=InterestingDataArray(jj,3); % get another object to compare
        if jj~=ii % if these are not same object, compare them
            if (abs(InterestingDataArray(jj,4)-InterestingDataArray(ii,4))) < 10 && (abs(InterestingDataArray(jj,5)-InterestingDataArray(ii,5)) < 10)
                % if these are the same object (within 10 pixels in each direction), place in vectors
                itn=itn+1;
                orig(itn)=currentobject; % original files
                dup(itn)=evalobj; % duplicates of original files
            end
        end
    end
end

% Go to airways (or Folder1) folder
pnchk=[pn '/objects/'];
duplicate_filenames = {};
if length(orig) >= 1
    for i = 1:length(orig)
       pairs{i} = sort([orig(i), dup(i)]);
    end

    if length(pairs) >= 1
        unique_pairs = {pairs{1}};
        i_ctr = 1;
        for i = 2:length(pairs)
            candidate = pairs{i};
            found = false;
            for j = 1:length(unique_pairs)
                if candidate == unique_pairs{j}
                    found = true;
                end
            end
            if found == false
                % Add to unique_pairs
                i_ctr = i_ctr + 1;
                unique_pairs{i_ctr} = candidate;
            end
        end
    end
    

    %% Open folder for checking duplicates
    % Check to see that the files exist

    for i = 1:length(unique_pairs)
        if exist([pnchk 'image_Object_' num2str(unique_pairs{i}(1)) '.mat'], 'file') == 0
            disp (['Help, the file associated to ' num2str(unique_pairs{i}(1)) ...
                " doesn't exist."])
        end
        if exist([pnchk 'image_Object_' num2str(unique_pairs{i}(2)) '.mat'], 'file') == 0
            disp (['Help, the file associated to ' num2str(unique_pairs{i}(2)) ...
                " doesn't exist."])
        end
        duplicate_filenames{end+1} = ['image_Object_' num2str(unique_pairs{i}(1)) '.mat'];
    end
end

if exist([pnchk 'Duplicates'])
    disp('The directory "Duplicates" exists');
else
    mkdir([pnchk 'Duplicates']);
end


savediary=1;
if savediary==1
    diary([pnchk 'Duplicates/test.log'])
    diary on
end


% Display paired duplicates
if length(orig)==0
    disp('No Duplicates');
else
    for kk=1:length(unique_pairs)
        disp(['Duplicates: ',...
            'Object ',num2str(unique_pairs{kk}(1)),' and ',...
            'object ',num2str(unique_pairs{kk}(2))]);
    end
end

% Turn off diary
savediary=1;
if savediary==1
    diary off
end


% move duplicates to "Duplicates" folder
for kk=1:length(duplicate_filenames)
    clear dupfn
    dupfn=duplicate_filenames{kk};
    if exist([pnchk dupfn])==2
        movefile([pnchk dupfn],[pnchk 'Duplicates/']);
    end
end

end