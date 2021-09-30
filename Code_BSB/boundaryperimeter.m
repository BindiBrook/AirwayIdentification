function perim = boundaryperimeter(coords,convrt)

perim=0;
coords=coords*convrt;
for i=1:size(coords,1)-1
    perim = perim+sqrt((coords(i+1,1)-coords(i,1))^2+(coords(i+1,2)-coords(i,2))^2);
end