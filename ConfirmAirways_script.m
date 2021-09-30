clear, 
close all, 
clc

% Manually identify airways out of all objects
dispmess=msgbox({'Now you need to check if the objects identified by ObjectIDscript.m are airways.'; ' '; 'You will be directed to select the tif you have analysed and then you will be shown each identified object in turn.'; ' ';'Please click OK to proceed'});
drawnow; waitfor(dispmess); clear dispmess
close;

cd Code_BSB
SelectAirways_FromMathsServers

%add in message please run: Analyse....
