function [level,mask]=Threshold_LAB(img,innerbound,outerbound,...
    innermask,outermask,transl,green)

global maskvar

defvals=[40 68];
level=defvals;

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

lowstep=.01; highstep=0.10;
lowminval=-127; lowmaxval=128;
highminval=0; highmaxval=128;

% GREY SLIDER
% OLD VERSION ent_grey = uicontrol('Style','pushbutton','Visible','on',...
%    'String','L.A.B. on "A" Channel',...
%    'Position',[ssz(3)/1.6 2*ssz(4)/3 105 30]);
ent_grey = uicontrol('Style','text','Visible','on',...
    'String','Threshold according to L.A.B. on "A" Channel. Use sliders to adjust upper and lower threshold, or use default values.',...
    'Position',[ssz(3)/1.6 2*ssz(4)/3 230 90]);
Gsld_low = uicontrol('Style', 'slider','Visible','on',...
    'Min',lowminval,'Max',lowmaxval,'Value',level(1),...
    'Sliderstep',[lowstep highstep],...
    'Position', [ssz(3)/1.85 ssz(4)-(ssz(4)/2)+50 200 20],.... % 90 20],...
    'Callback', @threshing_image,'Tag','Grey Low');
ed_Glow=uicontrol('Style','edit','Visible','on','String','0','Value',1,...
    'Position',[ssz(3)/1.85 ssz(4)-(ssz(4)/2)+100 90 20]);
Gsld_up = uicontrol('Style', 'slider','Visible','on',...
    'Min',highminval,'Max',highmaxval,'Value',level(2),...
    'Sliderstep',[lowstep highstep],...
    'Position', [3*ssz(3)/4.15 ssz(4)-(ssz(4)/2)+50 200 20],.... % 90 20],...
    'Callback', @threshing_image,'Tag','Grey Up'); % @adjustlevel);
ed_Gup=uicontrol('Style','edit','Visible','on','String','0','Value',1,...
    'Position',[3*ssz(3)/4.15 ssz(4)-(ssz(4)/2)+100 90 20]);

btn = uicontrol('Style', 'pushbutton','String','Apply Default Levels',...
        'Position', [ssz(3)/1.85 ssz(4)-(ssz(4)/2) 150 20],...
        'Callback', @applydefaults); 

btn = uicontrol('Style', 'pushbutton','String','OK',...
        'Position', [ssz(3)-(ssz(3)/4) ssz(4)-(ssz(4))+100 90 20],...
        'Callback', @onoff); 

uiwait

    %%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function onoff(object,~)
        valonoff=get(object,'value');
        if valonoff==1; close; end
    end
    
    function applydefaults(object,~)
        valsetdef=get(object,'value');
        if valsetdef==1
            level=defvals;
            set(ed_Glow,'String',level(1));
            set(ed_Gup,'String',level(2));
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
        % Map to L.A.B. space, threshold on the "A" channel
        % negative=imcomplement(rgb2gray(img));
        imgsz=size(img);
        LAB=rgb2lab(img);
        A=LAB(:,:,2);
        % Aimg(:,:,1)=zeros(imgsz(1),imgsz(2));
        % Aimg(:,:,2)=A;
        % Aimg(:,:,3)=zeros(imgsz(1),imgsz(2));
        
        out = A > Tvalue(1) & A < Tvalue(2); 
        
%         clear negative
%         dilated=bwmorph(out,'dilate',1);
%         clear out
%         de=bwmorph(dilated,'erode',1);
%         clear dilated
%         clean=bwareaopen(de,100);
%         clear de
%         fill=imfill(clean,'holes');
%         clear clean
        fill=out; clear out; % <- this is new, skip the above for now
        clear msk; msk=zeros(imgsz(1),imgsz(2));
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
    end
mask=maskvar; clear maskvar
end