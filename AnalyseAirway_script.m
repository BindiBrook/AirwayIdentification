clear, 
close all, 
clc,

% Draw BM and trace stain area
dispmess=msgbox({'Now, we need to determine the location of the basement membrane and relevant stain (aSMA or PSR) for all the objects you have identified as airways and are now in the "Airways" folder.';' '; 'Please click OK to proceed'});
drawnow; waitfor(dispmess); clear dispmess
close;
%add note about navigating to tif in the above message

%=======================================
[file,path]=uigetfile('*.tif','Select original tif file corresponding to the airways you want to analyse');
if path==0
    return
else
inputfile = [path file];
outputdir = [path file(1:end-4)];
end

cd Code_BSB
stain = histomenu('Is this an SMA or PSR stained image?','SMA','PSR');

if stain == 1
    SMAorPSR = 'SMA';
else
    SMAorPSR = 'PSR';
end
%===========================

Main_FromMathsServers(outputdir,file,SMAorPSR)
