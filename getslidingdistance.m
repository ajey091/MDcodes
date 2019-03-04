clear
close all
tic
c=linspecer(8);

fid = fopen('/Users/b119user/Downloads/Ajey/MDSimulations/GB9/tension/dump.xyz','r');
numatoms = cell2mat(textscan(fid,'%f\n',1,'HeaderLines',3));
fclose(fid);
fid = fopen('/Users/b119user/Downloads/Ajey/MDSimulations/GB9/tension/dump.xyz','r');
for ii = 1:1
    rawdata = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f\n',numatoms,'HeaderLines',9);
    data = cell2mat(rawdata);
    GBline = data(data(:,2)==4,:);
    p = polyfit(GBline(:,3),GBline(:,4),1);
    f1 = polyval(p,data(:,3));
    alllineind = mintersect(find(data(:,2)==3),find(data(:,5)>-10),find(data(:,5)<10));
    alllinedata = data(alllineind,:);
    indleft = find(p(1)*alllinedata(:,3) - alllinedata(:,4) + p(2) > 0);
    indright = find(p(1)*alllinedata(:,3) - alllinedata(:,4) + p(2) < 0);
    leftpoints = alllinedata(indleft,:);
    rightpoints = alllinedata(indright,:);
    R = [cosd(atand(p(1))) -sind(atand(p(1))); sind(atand(p(1))) cosd(atand(p(1)))];
    rotleftpoints = (R.'*leftpoints(:,3:4).').';
    rotleftpoints = rotleftpoints(intersect(find(rotleftpoints(:,2)<200),find(rotleftpoints(:,2)>-200)),:);
    line1indleft = intersect(find(rotleftpoints(:,1)<50), find(rotleftpoints(:,1)>-50));
    line2indleft = intersect(find(rotleftpoints(:,1)<-50), find(rotleftpoints(:,1)>-150));
    line3indleft = intersect(find(rotleftpoints(:,1)<150), find(rotleftpoints(:,1)>50));
    rotrightpoints = (R.'*rightpoints(:,3:4).').';
    rotrightpoints = rotrightpoints(intersect(find(rotrightpoints(:,2)<200),find(rotrightpoints(:,2)>-200)),:);
    line1indright = intersect(find(rotrightpoints(:,1)<50), find(rotrightpoints(:,1)>-50));
    line2indright = intersect(find(rotrightpoints(:,1)<-50), find(rotrightpoints(:,1)>-150));
    line3indright = intersect(find(rotrightpoints(:,1)<150), find(rotrightpoints(:,1)>50));    
    slidedist1(ii,1) = getperpdistance(rotleftpoints(line1indleft,:),rotrightpoints(line1indright,:));
    slidedist2(ii,1) = getperpdistance(rotleftpoints(line2indleft,:),rotrightpoints(line2indright,:));
    slidedist3(ii,1) = getperpdistance(rotleftpoints(line3indleft,:),rotrightpoints(line3indright,:));
    meanslidingdistwo(ii,1) = mean([slidedist1(ii,1),slidedist2(ii,1),slidedist3(ii,1)]);
end
fclose(fid);

fid = fopen('/Users/b119user/Downloads/Ajey/MDSimulations/GB9/tension_void2/dump.xyz','r');
numatoms = cell2mat(textscan(fid,'%f\n',1,'HeaderLines',3));
fclose(fid);
fid = fopen('/Users/b119user/Downloads/Ajey/MDSimulations/GB9/tension_void2/dump.xyz','r');
for ii = 1:300
    rawdata = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f\n',numatoms,'HeaderLines',9);
    data = cell2mat(rawdata);
    GBline = data(data(:,2)==4,:);
    p = polyfit(GBline(:,3),GBline(:,4),1);
    f1 = polyval(p,data(:,3));
    alllineind = mintersect(find(data(:,2)==3),find(data(:,5)>-10),find(data(:,5)<10));
    alllinedata = data(alllineind,:);
    indleft = find(p(1)*alllinedata(:,3) - alllinedata(:,4) + p(2) > 0);
    indright = find(p(1)*alllinedata(:,3) - alllinedata(:,4) + p(2) < 0);
    leftpoints = alllinedata(indleft,:);
    rightpoints = alllinedata(indright,:);
    R = [cosd(atand(p(1))) -sind(atand(p(1))); sind(atand(p(1))) cosd(atand(p(1)))];
    rotleftpoints = (R.'*leftpoints(:,3:4).').';
    rotleftpoints = rotleftpoints(intersect(find(rotleftpoints(:,2)<200),find(rotleftpoints(:,2)>-200)),:);
    line1indleft = intersect(find(rotleftpoints(:,1)<50), find(rotleftpoints(:,1)>-50));
    line2indleft = intersect(find(rotleftpoints(:,1)<-50), find(rotleftpoints(:,1)>-150));
    line3indleft = intersect(find(rotleftpoints(:,1)<150), find(rotleftpoints(:,1)>50));
    rotrightpoints = (R.'*rightpoints(:,3:4).').';
    rotrightpoints = rotrightpoints(intersect(find(rotrightpoints(:,2)<200),find(rotrightpoints(:,2)>-200)),:);
    line1indright = intersect(find(rotrightpoints(:,1)<50), find(rotrightpoints(:,1)>-50));
    line2indright = intersect(find(rotrightpoints(:,1)<-50), find(rotrightpoints(:,1)>-150));
    line3indright = intersect(find(rotrightpoints(:,1)<150), find(rotrightpoints(:,1)>50));    
    slidedist1(ii,1) = getperpdistance(rotleftpoints(line1indleft,:),rotrightpoints(line1indright,:));
    slidedist2(ii,1) = getperpdistance(rotleftpoints(line2indleft,:),rotrightpoints(line2indright,:));
    slidedist3(ii,1) = getperpdistance(rotleftpoints(line3indleft,:),rotrightpoints(line3indright,:));
    meanslidingdistw(ii,1) = mean([slidedist1(ii,1),slidedist2(ii,1),slidedist3(ii,1)]);
end
fclose(fid);

cmap = linspecer(2);
xq = linspace(0,300,20);
noslipdist = interp1(linspace(0,300,300),meanslidingdistwo,xq);
withslipdist = interp1(linspace(0,300,300),meanslidingdistw,xq);
noslipdist(noslipdist<0) = 0;
withslipdist(withslipdist<0) = 0;

figure('pos',[10,10,1000,800]),
plot(xq,noslipdist,'LineWidth',5,'color',cmap(1,:))
hold on
plot(xq,withslipdist,'LineWidth',5,'color',cmap(2,:))
xlabel('Macroscopic strain')
ylabel('Sliding Distance (A)')
legend('No slip interaction','With slip interaction','location','northwest')
set(gca,'FontSize',30)
save('GB9slidingdistances_new.mat','noslipdist','withslipdist')

% pnoslip = polyfit(xq,noslip,1)

toc
