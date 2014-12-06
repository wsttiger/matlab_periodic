clear all;
close all;

NL = 5;
maxdist = 15;

% lattice vectors
avec = [10.04579    0.00000    0.00000;
         0.00000   10.73176    0.00000;
         0.00000    0.00000   14.38271];
avec = avec';

% atoms
Y = [0.9793  0.0729  0.2500;
     0.0207  0.9271  0.7500;
     0.5207  0.5729  0.2500;
     0.4793  0.4271  0.7500];
O = [0.1210  0.4580  0.2500;
     0.8790  0.5420  0.7500;
     0.3790  0.9580  0.2500;
     0.6210  0.0420  0.7500;
     0.6910  0.3100  0.0580;
     0.3090  0.6900  0.9420;
     0.8090  0.8100  0.0580;
     0.1910  0.1900  0.9420;
     0.1910  0.1900  0.5580;
     0.8090  0.8100  0.4420;
     0.3090  0.6900  0.5580;
     0.6910  0.3100  0.4420];
 Ti = [0.0000  0.5000  0.0000;
       0.5000  0.0000  0.0000;
       0.5000  0.0000  0.5000;
       0.0000  0.5000  0.5000];

%atomslc = [Y; O; Ti];
atomslc = [Ti];
% lattice coordinates to cartesian coordinates
atomsrc = (avec*atomslc')';

% generate lots of atoms
atomssc_t1 = zeros((NL+1)*(NL+1)*(NL+1),3);
tindx = 1;
for i1 = -NL:NL
    for i2 = -NL:NL
        for i3 = -NL:NL
            for ias = 1:size(atomslc,1)
                atomssc_t1(tindx,:) = atomslc(ias,:) + [i1 i2 i3];
                tindx = tindx + 1;
            end
        end
    end
end
% lattice coordinates to cartesian coordinates
atomssc_t2 = (avec*atomssc_t1')';
% sort by the norm
atomssc_t2_norm = zeros(size(atomssc_t2,1),1);
for i = 1:size(atomssc_t2,1)
    atomssc_t2_norm(i) = norm(atomssc_t2(i,:));
end
atomssc_t3 = sortrows([atomssc_t2 atomssc_t2_norm],4);
% filter by maxdist
maxi = find(atomssc_t3(:,4) > maxdist);
atomssc = atomssc_t3(1:maxi,:);
