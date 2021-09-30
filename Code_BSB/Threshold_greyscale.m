function [f,level,mask]=Threshold_greyscale(img,innerbound,outerbound,...
    innermask,outermask,transl,green)

global T1 maskvar

level=[100 120];

ssz = get(0,'ScreenSize');

% [left bottom width height]
hfig=figure('Visible','on','Name','IMAGE THRESHOLDING',...
    'NumberTitle','off',...
    'Position',[ssz(3)/90 ssz(4)/90 ssz(3)-(ssz(3)/8) ssz(4)-(ssz(4)/6)]);
hAxes=axes('Parent',hfig,...
    'Position',[-0.1 0.1 0.75 0.75]); 

imshow(img,'Parent',hAxes); 
    hold on;
    plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
    hold on;
    plot(outerbound(:,2),outerbound(:,1),'b','LineWidth',3);

range=.01; minval=0; maxval=256;

% GREY SLIDER
%SAME EDIT AS IN ..._LAB. SIZING ETC NEEDS TESTING
%ent_grey = uicontrol('Style','pushbutton','Visible','on',...
%    'String','GREYSCALE',...
%    'Position',[ssz(3)-(ssz(3)/2.75) ssz(4)-(ssz(4)/2)+300 105 30]);
ent_grey = uicontrol('Style','pushbutton','Visible','on',...
    'String','Threshold according to greyscale image data. Use sliders to adjust upper and lower threshold, or use default values.',...
    'Position',[ssz(3)-(ssz(3)/2.75) ssz(4)-(ssz(4)/2)+300 230 90]);
Gsld_low = uicontrol('Style', 'slider','Visible','on',...
    'Min',minval,'Max',maxval,'Value',level(1),'Sliderstep',[range range],...
    'Position', [ssz(3)-(ssz(3)/2.75) ssz(4)-(ssz(4)/2)+200 200 20],.... % 90 20],...
    'Callback', @threshing_image,'Tag','Grey Low');
ed_Glow=uicontrol('Style','edit','Visible','on','String','0','Value',1,...
    'Position',[ssz(3)-(ssz(3)/2.75) ssz(4)-(ssz(4)/2)+250 90 20]);
Gsld_up = uicontrol('Style', 'slider','Visible','on',...
    'Min',minval,'Max',maxval,'Value',level(2),'Sliderstep',[range range],...
    'Position', [ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4)/2)+200 200 20],.... % 90 20],...
    'Callback', @threshing_image,'Tag','Grey Up'); % @adjustlevel);
ed_Gup=uicontrol('Style','edit','Visible','on','String','0','Value',1,...
    'Position',[ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4)/2)+250 90 20]);

btn = uicontrol('Style', 'pushbutton','String','OK',...
        'Position', [ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4))+250 90 20],...
        'Callback', @onoff); 
    
uiwait

    %%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function onoff(object,~)
        valonoff=get(object,'value');
        if valonoff==1; close; end
    end
    
    function threshing_image(object,~)
        if strcmpi(get(object, 'Tag'), 'Grey Low')
            level(1)=get(object,'value');
            set(ed_Glow,'String',level(1));
        elseif strcmpi(get(object, 'Tag'), 'Grey Up')
            level(2)=get(object,'value');
            set(ed_Gup,'String',level(2));
        end
        maskvar=Imthreshold_PSR(level);
        imshow(img);
            hold on;
            plot(innerbound(:,2),innerbound(:,1),'r','LineWidth',3);
            hold on;
            plot(outerbound(:,2),outerbound(:,1),'b','LineWidth',3);
            hold on;
            h = imshow(green);
            hold off;
            set(h,'AlphaData',maskvar);
    end

    function msk=Imthreshold_PSR(Tvalue)
        %Replace the pixel value either with 0 or 255
        % Convert to grayscale, invert, normalise to [0 255]
        negative=imcomplement(rgb2gray(img));
        
        out = negative > Tvalue(1) & negative < Tvalue(2); 
        
        clear negative
        dilated=bwmorph(out,'dilate',1);
        clear out
        de=bwmorph(dilated,'erode',1);
        clear dilated
        clean=bwareaopen(de,100);
        clear de
        fill=imfill(clean,'holes');
        clear clean
        clear msk; imgsz=size(img); msk=zeros(imgsz(1),imgsz(2));
        for i=1:imgsz(1)
            for j=1:imgsz(2)
                if fill(i,j)==1
                    msk(i,j)=transl;
                else
                end
            end
        end
        msk=msk.*(double(innermask==0));
        msk=msk.*(double(outermask==1));
        
        % Im=img;
        % Im(:,:,2)=uint16(double(img(:,:,2)).*(1-maskvar)+...
        %     255*maskvar); % make green overlay
    end
        
f=T1; clear T1
mask=maskvar; clear maskvar
end