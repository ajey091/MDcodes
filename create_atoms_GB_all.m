clear
close all
tic

% This code creates a bicrystal given the euler angles of the two grains making up the bicrystal, and the grain boundary normal vector
% Then, it iterates through the several atom configurations at the GB and identifies the lowest energy structure by minimizaion runs on these structures.


latconst = 4.04;

boxsize_x = 400;
boxsize_y = 400;
boxsize_z = boxsize_x;
GBnormvect = [0.646,0.7327,0]; % GB normal vector
GBnormvect = GBnormvect./norm(GBnormvect);
GBtangvect = cross(GBnormvect,[0,0,1]);

rotmat1 = [360,97,93]*pi/180; % Euler angles of orientation of grain A
rotmat2 = [358,70,96]*pi/180; % Euler angles of orientation of grain B
csrotmat1 = bungeRotationSample2Crystal(rotmat1).';
csrotmat2 = bungeRotationSample2Crystal(rotmat2).';

[crystalposx,crystalposy,crystalposz] = meshgrid([-boxsize_x:latconst:boxsize_x,-boxsize_x:latconst:boxsize_y,-boxsize_x:latconst:boxsize_z]);

crystalposallunique1 = [crystalposx(:),crystalposy(:),crystalposz(:)];
crystalposallunique2 = [crystalposx(:)+0.5*latconst,crystalposy(:)+0.5*latconst,crystalposz(:)];
crystalposallunique3 = [crystalposx(:),crystalposy(:)+0.5*latconst,crystalposz(:)+0.5*latconst];
crystalposallunique4 = [crystalposx(:)+0.5*latconst,crystalposy(:),crystalposz(:)+0.5*latconst];

crystalposallunique = [unique(crystalposallunique1,'rows');
    unique(crystalposallunique2,'rows');
    unique(crystalposallunique3,'rows');
    unique(crystalposallunique4,'rows')];

crystal1samp = (csrotmat1*crystalposallunique.').';
crystal2samp = (csrotmat2*crystalposallunique.').';

midcrys1 = mean(crystal1samp,1);
midcrys2 = mean(crystal2samp,1);
GBpoint = mean([midcrys1;midcrys2]);
syms x y z

GBplane_eq = GBnormvect(1)*(x-GBpoint(1)) + GBnormvect(2)*(y-GBpoint(2)) + GBnormvect(3)*(z-GBpoint(3));

[A,B] = equationsToMatrix(GBplane_eq,[x,y,z]);
A = eval(A);
B = eval(B);

vect1 = [A,B];
vect2 = [crystal1samp,ones(length(crystal1samp),1)];
planeval1 = sum(vect1.'.*vect2.');
planeval1ind = find(planeval1 > 1);
crystal1samp_snipped = crystal1samp(planeval1ind,:);

vect1 = [A,B];
vect2 = [crystal2samp,ones(length(crystal2samp),1)];
planeval2 = sum(vect1.'.*vect2.');
planeval2ind = find(planeval2 < -1);
% 
trans1 = linspace(0,latconst*2,20);
trans2 = linspace(0,latconst*2,20);

%%
count = 0;
figure('pos',[10,10,1000,600]),

for ii = 1:2
    for jj = 1:2
        count = count + 1;
        crystal2samp_snipped = crystal2samp(planeval2ind,:) + GBtangvect.*trans1(ii) + [0,0,1]*trans2(jj);
        
        finalboxleftlim = -200;
        finalboxrightlim = 200;
	finalboxtoplim =    200;
	finalboxbottomlim = -200;
	finalboxthick1 = 20;
	finalboxthick2 = -20;
        
        finalinds1 = mintersect(find(crystal1samp_snipped(:,1) > finalboxleftlim), ...
            find(crystal1samp_snipped(:,1) < finalboxrightlim),...
            find(crystal1samp_snipped(:,2) < finalboxtoplim),...
            find(crystal1samp_snipped(:,2) > finalboxbottomlim),...
            find(crystal1samp_snipped(:,3) < finalboxthick1),...
            find(crystal1samp_snipped(:,3) > finalboxthick2));
        finalinds2 = mintersect(find(crystal2samp_snipped(:,1) > finalboxleftlim), ...
            find(crystal2samp_snipped(:,1) < finalboxrightlim),...
            find(crystal2samp_snipped(:,2) < finalboxtoplim),...
            find(crystal2samp_snipped(:,2) > finalboxbottomlim),...
            find(crystal2samp_snipped(:,3) < finalboxthick1),...
            find(crystal2samp_snipped(:,3) > finalboxthick2));
        
        filename = sprintf('datafile_GB33_%d_%d',ii,jj);
        fid = fopen(filename,'w');
        
        fprintf(fid,'Position data for GB33\n\n');
        fprintf(fid,'%d atoms\n',length(finalinds1)+length(finalinds2));
        fprintf(fid,'2 atom types\n');
        fprintf(fid,'%f %f xlo xhi\n',finalboxleftlim-1,finalboxrightlim+1);
        fprintf(fid,'%f %f ylo yhi\n',finalboxbottomlim-1,finalboxtoplim+1);
        fprintf(fid,'%f %f zlo zhi\n\n',finalboxthick2-1,finalboxthick1+1);
        fprintf(fid,'0.0 0.0 0.0 xy xz yz\n\n');
        fprintf(fid,'Atoms\n\n');
        
        fprintf(fid,'%d %d %f %f %f\n',[[1:length(finalinds1)].',ones(length(finalinds1),1),crystal1samp_snipped(finalinds1,1),...
            crystal1samp_snipped(finalinds1,2),crystal1samp_snipped(finalinds1,3)].');
        
        fprintf(fid,'%d %d %f %f %f\n',[[length(finalinds1)+1:length(finalinds1)+length(finalinds2)].',2.*ones(length(finalinds2),1),crystal2samp_snipped(finalinds2,1),...
            crystal2samp_snipped(finalinds2,2),crystal2samp_snipped(finalinds2,3)].');
        
        fclose(fid);
        system(sprintf('sed  "s/datafile_GB33/datafile_GB33_%d_%d/g" Al_GB33.in > Al_GB33_temp.in',ii,jj));
        
	system(sprintf('mpirun -np 80 lmp_mpi < Al_GB33_temp.in > log.lammps_%d_%d',ii,jj));
    	[status,cmdout] = system(sprintf('grep -i "Loop time" -B 1 log.lammps_%d_%d',ii,jj));
        poten(count) = str2double(cmdout(12:22));
        [~,cmdout2] = system(sprintf('tail -n 2 log.lammps_%d_%d | head -1',ii,jj));
        atoms(count) = str2double(cmdout2);
        Eatom(count) = poten(count)/atoms(count);

     end
end
plot(1:length(Eatom(:)),Eatom(:))
hold on

xlabel('Input structure ID')
xlim([1,400])
ylabel('Converged potential energy per atom (eV)')
set(gca,'FontSize',26)

toc
