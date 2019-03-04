function val = getperpdistance (coords1,coords2)
    mid1 = [mean(coords1(:,1)),mean(coords1(:,2))];
    mid2 = [mean(coords2(:,1)),mean(coords2(:,2))];
    p1 = polyfit(coords1(:,1),coords1(:,2),1);
    p2 = polyfit(coords2(:,1),coords2(:,2),1);
    xvals = linspace(-200,200,100);
    perpline1 = -1/p1(1) * (xvals - mid1(1)) + mid1(2);
    parallelline1 = p1(1) * (xvals - mid1(1)) + mid1(2); 
    parallelline2 = p2(1) * (xvals - mid2(1)) + mid2(2); 
    [x1,y1] = intersections(xvals,perpline1,xvals,parallelline1);
    [x2,y2] = intersections(xvals,perpline1,xvals,parallelline2);
    X = [x1,y1;x2,y2];
    val = pdist(X);
%     figure,
%     plot(coords1(:,1),coords1(:,2),'.')
%     hold on
%     plot(coords2(:,1),coords2(:,2),'.')
%     plot(mid1(:,1),mid1(:,2),'.')
%     plot(xvals,perpline1)
%     plot(xvals,parallelline2)
%     plot(xvals,parallelline1)
%     axis equal
end
