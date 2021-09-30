function outerbound = boundarylength(innermasknew)

% Find points along boundaries
j=1; %initialize j
i=1; %initialize i
for i=1:length(innermasknew(1,:))
    m=0; %initialize m
    for j=1:length(innermasknew(:,1))
        m=m+1;
        lencol(i)=m;
        if innermasknew(j,i)~=0, break, end
    end
end
j=min(lencol);
i=find(lencol==j);
i=i(1); % Only use one value for starting position, i

% Trace boundary
outerbound=bwtraceboundary(innermasknew,[j i],'N');
% plot(outerbound(:,2),outerbound(:,1),'g','LineWidth',3);

end