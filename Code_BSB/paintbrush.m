function [tempmask]=paintbrush(A,tempmask,transl,green)

ssz = get(0,'ScreenSize');
% [left bottom width height]
figure('Visible','on','Name','MANUALLY SELECT STAIN AREAS',...
    'NumberTitle','off',...
    'Position',[ssz(3)/50 ssz(4)/50 ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4)/3)]);
% axes('units','pixels',...
%     'Position',[ssz(3)/45 ssz(4)/45 ssz(3)-(ssz(3)/5) ssz(4)-(ssz(4)/4)]); % [ssz(3)/45 ssz(4)/45 ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4)/3)]);
hold on;

sz=size(A); 
greenmask(:,:,1)=zeros(sz(1),sz(2));
greenmask(:,:,2)=255*tempmask;
greenmask(:,:,3)=zeros(sz(1),sz(2));
applymask(:,:,1)=ones(sz(1),sz(2)).*imcomplement(tempmask);
applymask(:,:,2)=ones(sz(1),sz(2)).*imcomplement(tempmask);
applymask(:,:,3)=ones(sz(1),sz(2)).*imcomplement(tempmask);
applymask=uint8(applymask);
% displayimg=(A.*applymask+uint8(greenmask)); 
imagesc(A); % imagesc(displayimg);
hold on;
h=imshow(green);
hold off;

set(h,'AlphaData',transl*greenmask(:,:,2)/255);

colormap('default'); % This sets the colormap
% colormap(gray);

% Button to export data to *.mat file
%THIS DOESN'T SEEM TO DO ANYTHING THE SAVE COMMAND IN THE BELOW FUNCTION IS
%COMMENTED OUT.
%h1 = uicontrol('Style', 'pushbutton', ...
%    'Position', [100 450 200 50], ...
%    'String','Save mask value',...
%    'Callback', @savemask);

% Slider to change brush size, w/ Tag "BrushSize"
h2 = uicontrol('Style', 'slider', ...
    'Position', [100 500 200 50], ...
    'Value',5,...
    'Max',30,...
    'Min',1,...
    'Tag','BrushSize',...
    'SliderStep',[0.2  0.2]);

h3 = uicontrol('Style', 'text', ...
     'Position', [100 550 200 50], ...
     'String',...
     ['Brush Size (slide to change brush size, left click on image to paint, ',...
     'right click to erase...'],...
     'Tag','txtErase');

% h4 = uicontrol('Style', 'togglebutton',...
%     'Position', [100 400 200 50], ...
%     'String','Push To Change to Brush Shape',...
%     'Min', 1,...
%     'Max', 0,...
%     'Tag','BrushShape');

% h5 = uicontrol('Style', 'text', ...
%      'Position', [100 450 200 50], ...
%      'String',...
%      [{'This button changes the brush shape: '},...
%      {'OFF (Light Grey) is Circle '},...
%      {'ON (Dark Grey) is Square'}],...
%      'Tag','txtErase');

btn = uicontrol('Style', 'pushbutton','String','OK',...
        'Position', [100 300 200 50],...
        'Callback', @onoff); 
    
axis image;

% Get size of the image
m=size(A,1);n=size(A,2);

% Set limits on figure??
xlim([1 n]);
ylim([1 m]);
 
% Unpack gui object
gui = get(gcf,'UserData'); % obtains user data from gui

% set(gcf,'PointerShapeCData',cc);
set(gcf,'Pointer','crosshair');

% Make a fresh figure window
set(gcf,'WindowButtonDownFcn',@startmovit);

% Store gui object
set(gcf,'UserData',{gui;A;applymask;tempmask});

uiwait

function onoff(object,~)
        val=get(object,'value');
        if val==1
            close;
        end
end

function savemask(src,evnt)
% Unpack gui object
temp = get(gcf,'UserData'); % This gets data from gui
gui=temp{1}; % This is data from gui wrapper
             % this function does not appear to use it
A=temp{2}; % This is the new image w/ the drawn/deleted pixels
applymask=temp{3};
tempmask=temp{4};
% save ('A.mat','A'); % This saves the final image to file
% save ('mask.mat','mask'); % This saves mask to file
end

function display(src,evnt)
% It looks like this function iterates through the slices
% Unpack gui object
temp = get(gcf,'UserData'); % This gets data from gui
gui=temp{1}; % This is data from gui wrapper
A=temp{2}; % This is the new image w/ the drawn/deleted pixels
applymask=temp{3};
tempmask=temp{4};
sz=size(A);
greenmask(:,:,1)=zeros(sz(1),sz(2));
greenmask(:,:,2)=255*tempmask;
greenmask(:,:,3)=zeros(sz(1),sz(2));
% displayimg=(A.*applymask+uint8(greenmask)); 
imagesc(A); % imagesc(displayimg);
hold on;
h = imshow(green);
hold off;
set(h,'AlphaData',transl*greenmask(:,:,2)/255);
colormap('default');
axis image;
end

function startmovit(src,evnt)
% disp('start tracing')
% It looks like this function draws/erases
temp = get(gcf,'UserData'); % This gets data from gui
gui=temp{1}; % This is data from gui wrapper
A=temp{2}; % This gets the image from the gui
applymask=temp{3}; % This gets the mask from gui
tempmask=temp{4};

% Get brush size from slider
brslider=findobj(gcf,'Style','Slider','-and','Tag','BrushSize'); % Brush size slider
brushsize=get(brslider,'Value'); % Brush size
% disp(num2str(brushsize))

pos = get(gca,'CurrentPoint'); % get current mouse position
rightleftclick=get(src,'SelectionType'); % Right or left mouse click
% Normal = left, alt=right

[m,n]=size(A); % get size of A, which I think is the image displacyed

% Get mouse position in integers
posm=round(pos(1,2));
posn=round(pos(1,1));

% Set draw/erase position
xrange=max(posm-brushsize/2,1):min(posm+brushsize/2,m);
yrange=max(posn-brushsize/2,1):min(posn+brushsize/2,n);

% Get info on brush shape
% togbutton=findobj(gcf,'Style','togglebutton','-and','Tag','BrushShape'); % Brush shape toggle
brushshape=1; % keep a circle % get(togbutton,'Value'); % Brush shape
% disp(num2str(brushshape))

if brushshape==1
    % Brush is a circle
    for i=1:length(xrange)
        for j=1:length(yrange)
            if sqrt((xrange(i)-posm)^2+(yrange(j)-posn)^2)<=brushsize/2
                if strcmp(rightleftclick,'normal') % Draw if left click
                    applymask(uint32(xrange(i)),uint32(yrange(j)),:)=0;
                    tempmask(uint32(xrange(i)),uint32(yrange(j)))=1;
                elseif strcmp(rightleftclick,'alt') % Erase if right click
                    applymask(uint32(xrange(i)),uint32(yrange(j)),:)=1;
                    tempmask(uint32(xrange(i)),uint32(yrange(j)))=0;
                end
            end
        end
    end
else
    % Brush is a square
    xrange=uint32(xrange);
    yrange=uint32(yrange);
    % Switch between draw and erase
    if strcmp(rightleftclick,'normal') % Draw if left click
        applymask(xrange,yrange,:)=0;
        tempmask(xrange,yrange)=1;
    elseif strcmp(rightleftclick,'alt') % Erase if right click
        applymask(xrange,yrange,:)=1;
        tempmask(xrange,yrange)=0;
    end
end

% Display in the window
sz=size(A);
greenmask(:,:,1)=zeros(sz(1),sz(2));
greenmask(:,:,2)=255*tempmask;
greenmask(:,:,3)=zeros(sz(1),sz(2));
% displayimg=(A.*applymask+uint8(greenmask)); 
imagesc(A); % imagesc(displayimg);
hold on;
h = imshow(green);
hold off;
set(h,'AlphaData',transl*greenmask(:,:,2)/255);
colormap('default');
axis image;

% Set callbacks
gui.currenthandle = src;
thisfig = gcbf();
set(thisfig,'WindowButtonUpFcn',@stopmovit);
 
% Place image and data back in the gui
set(gcf,'UserData',{gui;A;applymask;tempmask});
end

function stopmovit(src,evnt)
% disp('stop tracing')
% Clean up the evidence ...
thisfig = gcbf();
ud = get(gcf,'UserData'); 
set(thisfig,'WindowButtonUpFcn','');
set(thisfig,'WindowButtonMotionFcn','');
axis image;
set(gcf,'UserData',ud);
end
end