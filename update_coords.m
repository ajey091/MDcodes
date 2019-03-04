clear
close all

tic

GBnormvect = [0.646,0.732,0];
GBnormvect = GBnormvect./norm(GBnormvect);

[k1,k2] = system('grep -n TIMESTEP GB33.xyz | tail -1 | cut -d: -f1');
linenum = str2num(k2);

fid = fopen('GB33.xyz','r');
coordlimits = textscan(fid,'%f %f %f\n',3,'HeaderLines',linenum+4);
fclose(fid);
coordlimits = cell2mat(coordlimits);

fid = fopen('GB33.xyz','r');
alldata = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f\n','HeaderLines',linenum+8);
fclose(fid);

data = cell2mat(alldata);

finalboxleftlim = -10*20;
finalboxrightlim = 10*20;

GBpoint = [ -0.9834   -0.8770    0.9596];

scribelines = linspace(finalboxleftlim,finalboxrightlim,5);

datalabel = ones(length(data),1);

for ii = 1:length(data)
    if (GBnormvect(1)*data(ii,3) + GBnormvect(2)*data(ii,4) > 2)
        datalabel(ii,1) = 2;
    end
end
for jj = 1:length(scribelines)
    for ii = 1:length(data)
        if (-GBnormvect(1)*data(ii,3) + GBnormvect(2)*data(ii,4) > scribelines(jj)-2 && -GBnormvect(1)*data(ii,3) + GBnormvect(2)*data(ii,4) < scribelines(jj)+2)
            datalabel(ii,1) = 3;
        end
        if (GBnormvect(1)*data(ii,3) + GBnormvect(2)*data(ii,4) > GBpoint(1)-4 && GBnormvect(1)*data(ii,3) + GBnormvect(2)*data(ii,4) < GBpoint(2)+4)
                datalabel(ii,1) = 4;
        end
    end 
end

num2 = 100;
num3 = 400;
num = 8;
voididx = mintersect(find(num*data(:,3)+data(:,4)+num3<-num2),find(num*data(:,3)-data(:,4)+num3>num2));
data(voididx,:) = [];
datalabel(voididx,:) = [];
% 
voididx = mintersect(find(num*data(:,3)+data(:,4)-num3>num2),find(num*data(:,3)-data(:,4)-num3<-num2));
data(voididx,:) = [];
datalabel(voididx,:) = [];


% Writing lammps input file
filename = sprintf('datafile_GB33_strain_minimized_4types_void.in');
fid = fopen(filename,'w');

fprintf(fid,'Position data for GB 33\n\n');
fprintf(fid,'%d atoms\n',length(data));
fprintf(fid,'4 atom types\n');
fprintf(fid,'%f %f xlo xhi\n',coordlimits(1,1)-1,coordlimits(1,2)+1);
fprintf(fid,'%f %f ylo yhi\n',coordlimits(2,1)-1,coordlimits(2,2)+1);
fprintf(fid,'%f %f zlo zhi\n\n',coordlimits(3,1)-1,coordlimits(3,2)+1);
fprintf(fid,'0.0 0.0 0.0 xy xz yz\n\n');
fprintf(fid,'Atoms\n\n');

fprintf(fid,'%d %d %f %f %f\n',[[1:length(data(datalabel==1))].',ones(length(data(datalabel==1)),1),...
    data(datalabel==1,3),data(datalabel==1,4),data(datalabel==1,5)].');

fprintf(fid,'%d %d %f %f %f\n',[[length(data(datalabel==1))+1:length(data(datalabel==1))+...
    length(data(datalabel==2))].',ones(length(data(datalabel==2)),1)*2,data(datalabel==2,3),...
    data(datalabel==2,4),data(datalabel==2,5)].');

fprintf(fid,'%d %d %f %f %f\n',[[length(data(datalabel==1))+length(data(datalabel==2))+...
    1:length(data(datalabel==1))+length(data(datalabel==2))+length(data(datalabel==3))].',...
    ones(length(data(datalabel==3)),1)*3,data(datalabel==3,3),data(datalabel==3,4),data(datalabel==3,5)].');

fprintf(fid,'%d %d %f %f %f\n',[[length(data(datalabel==1))+length(data(datalabel==2))+...
    length(data(datalabel==3))+1:length(data(datalabel==1))+length(data(datalabel==2))+...
    length(data(datalabel==3))+length(data(datalabel==4))].',ones(length(data(datalabel==4)),1)*4,...
    data(datalabel==4,3),data(datalabel==4,4),data(datalabel==4,5)].');


fclose(fid);


toc
