[file,path]=uigetfile('*.tif','Select the tif file you want to analyse');
if path==0
    return
else
inputfile = [path file];
outputdir = [path file(1:end-4)];

cd Code_BSB
stain = histomenu('Is this an SMA or PSR stained image?','SMA','PSR');

if stain == 1
    SMAorPSR = 'SMA';
else
    SMAorPSR = 'PSR';
end

% Analyse tiff

dispmess=msgbox('The tiff will be analysed next. This might take a few hours. Please click OK to proceed.');
drawnow; waitfor(dispmess); clear dispmess
close;

ObjectID_loop(inputfile,'output_directory',outputdir,'SMAorPSR',SMAorPSR,'tile_size',10000)

dispmess=msgbox('All done! Next, you will need to confirm that the found objects are (or are not) airways. Please run "ConfirmAirways_script.m". Please click OK to proceed.');
drawnow; waitfor(dispmess); clear dispmess
close;

end
cd ../