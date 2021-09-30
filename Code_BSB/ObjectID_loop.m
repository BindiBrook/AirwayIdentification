function ObjectID_loop(varargin)

    p = inputParser;
    % the first argument is required, the filename of the image
    addRequired(p, 'input_filename', @(x)validateattributes(x,{'char'},{'nonempty'}))
    % output_directory defaults to pwd
    cd ../
    addParameter(p, 'output_directory', pwd, @(x)validateattributes(x,{'char'},{'nonempty'}))
    cd Code_BSB
    % tile_size has default 5000
    addParameter(p, 'tile_size', 5000, @(x)validateattributes(x,{'numeric'},{'nonnegative'}))
    % overlap_size has default -1 here to indicate that no value has been 
    % passed. In that case, the default (based upon PixelsToMicrons) will 
    % be calcualated and set inside Preprocessor_MathsServers_Func().
    addParameter(p, 'overlap_size', -1, @(x)validateattributes(x,{'numeric'},{'nonnegative'}))
    % SMAorPSR can take value 'SMA' or 'PSR', and defaults to 'SMA'
    addParameter(p, 'SMAorPSR', 'SMA', @(x)validateattributes(x,{'char'},{'scalartext'}))
    parse(p,varargin{:})

    % If tile size is supplied as a scalar, then use a tile of that size
    % square.
    if size(p.Results.tile_size, 2) == 2
        tile_size = p.Results.tile_size;
    else
        tile_size = [p.Results.tile_size, p.Results.tile_size];
    end
    % If overlap size is supplied as a scalar, then use an overlap of that size
    % in both dimensions.
    if size(p.Results.overlap_size, 2) == 2
        overlap_size = p.Results.overlap_size;
    else
        overlap_size = [p.Results.overlap_size, p.Results.overlap_size];
    end

% GET CURRENT FILE/PATH
[~, file, ~] = fileparts(p.Results.input_filename);

    % Create subdirectories for saving results if does not exist already
    full_output_dir = [p.Results.output_directory];
  
    if isfolder(full_output_dir)
        rmdir(full_output_dir,'s')
        mkdir(full_output_dir)
    else
        mkdir(full_output_dir)
    end
    
    % Call program
     returnVal = Preprocessor_MathsServers_Func(p.Results.input_filename,full_output_dir,p.Results.SMAorPSR,tile_size,overlap_size);
     if returnVal == 1
         return
     end
  
    % Eliminate duplicates
    CheckForDuplicates(full_output_dir, 'InfoFile.mat');
    
    % Filter objects
    if strcmp(p.Results.SMAorPSR,'SMA') % Image is of SMA
        SMAPSRindex=1;
    else
        SMAPSRindex=2; % Image is of PSR
    end
    FilterObjectFiles(full_output_dir,SMAPSRindex);
end
