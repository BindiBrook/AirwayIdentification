function [bzbound,sz,loc]=bezier(I,boundary,halfwid)

figure(1)
imshow(I)
hold on
plot(boundary(:,2),boundary(:,1),'g','LineWidth',3)
title('Mouse-click or touch-screen to change boundary position','FontWeight','normal')
[x,y]=ginput(1);

sz=size(boundary);
%disp(['length of boundary is ',num2str(sz(1))])
ptmid=[boundary(1,2) boundary(1,1)];
for k=1:sz(1)
    if sqrt((boundary(k,1)-y)^2+(boundary(k,2)-x)^2) < sqrt((ptmid(2)-y)^2+(ptmid(1)-x)^2)
        ptmid=[boundary(k,2) boundary(k,1)];
        loc=k;
    end
end
%disp(['loc is ',num2str(loc)])

%figure; imshow(I)
hold on;
plot(boundary(:,2),boundary(:,1),'g',...,
    ptmid(1),ptmid(2),'r.','LineWidth',3,'MarkerSize',10);

% halfwid=10; % pixels
if 2*halfwid>sz(1)
    halfwid=floor(sz(1)/2);
    disp(['Error!!! Chosen smoothness parameter is greater than size of boundary, ',...
        're-setting to half size of boundary, ',...
        'smoothness = ',num2str(halfwid)]);
end

flag=0;
if loc-halfwid<1 || loc+halfwid>sz(1)
    flag=1;
    % Rotate the boundary by 180 degrees
    halfpt=floor(sz(1)/2);
    for k=1:sz(1)
        if k+halfpt>sz(1)
            rotbound(k,1)=boundary(k+halfpt-sz(1),1);
            rotbound(k,2)=boundary(k+halfpt-sz(1),2);
        else
            rotbound(k,1)=boundary(k+halfpt,1);
            rotbound(k,2)=boundary(k+halfpt,2);
        end
    end
    clear boundary; boundary=rotbound;
    if loc-halfwid<1
        pos1=sz(1)+(loc-halfwid);
    else
        pos1=loc-halfwid;
    end
    if loc+halfwid>sz(1)
        pos2=(loc+halfwid)-sz(1);
    else
        pos2=loc+halfwid;
    end
    pos1=pos1-halfpt;
    pos2=pos2+halfpt;
else
    pos1=loc-halfwid;
    pos2=loc+halfwid;
end

% if loc-halfwid<1
%     pos1=sz(1)+(loc-halfwid);
% else
%     pos1=loc-halfwid;
% end
% if loc+halfwid>sz(1);
%     pos2=(loc+halfwid)-sz(1);
% else
%     pos2=loc+halfwid;
% end
%disp(['pos1 is ',num2str(pos1)])
%disp(['pos2 is ',num2str(pos2)])

pt1=[boundary(pos1,2) boundary(pos1,1)]';
pt3=[boundary(pos2,2) boundary(pos2,1)]';
pt2=[x y]';
% disp(['loc - halfwid = ',num2str(loc-halfwid)])

% figure; 
% imshow(I)
% hold on;
% plot(boundary(:,2),boundary(:,1),'g',...,
%     ptmid(1),ptmid(2),'r.',...
%     pt1(1),pt1(2),'b.',...
%     pt2(1),pt2(2),'c.',...
%     pt3(1),pt3(2),'r.',...
%     'LineWidth',3,'MarkerSize',20);

t=linspace(0,1,101);
pts = kron((1-t).^2,pt1) + kron(2*(1-t).*t,pt2) + kron(t.^2,pt3);

% hold on
% plot(pts(1,:),pts(2,:))
% hold off

% Now stitch bezier curve back into boundary
sw=0;
it=1;
% newbound=zeros(2,sz(1)+length(pts(1,:)));
for k=1:sz(1)
    if boundary(k,2)==pt1(1) && boundary(k,1)==pt1(2)
        sw=1;
        for kk=1:length(pts(1,:))
            bzbound(it,2)=pts(1,kk); bzbound(it,1)=pts(2,kk);
            it=it+1;
        end
    end
    if boundary(k,2)==pt3(1) && boundary(k,1)==pt3(2)
        sw=0;
    end
    if sw==0
        bzbound(it,2)=boundary(k,2); bzbound(it,1)=boundary(k,1);
        it=it+1;
    end
end

%% WORK ON THIS BIT!!
if flag==1
    clear rotbound
    top=find(bzbound(:,2)==max(bzbound(:,2)));
    it=0;
    for k=top(1):length(bzbound(:,1))
        it=it+1;
        rotbound(it,1)=bzbound(k,1);
        rotbound(it,2)=bzbound(k,2);
    end
    for k=1:top(1)-1
        it=it+1;
        rotbound(it,1)=bzbound(k,1);
        rotbound(it,2)=bzbound(k,2);
    end
    clear bzbound; bzbound=rotbound;
end
end