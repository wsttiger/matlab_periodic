clear all;
close all;

thop = 0.5;

% original lattice vectors
avec = [1.0 0.0;
        0.0 1.0];


% atoms in the original unit cell
atomslc = [0.0, 0.0];

% generate lots of atoms (also record the lattice vectors to get
% these atoms)
NL = 30;
tindx = 1;
atomssc_t1 = zeros((NL+1)*(NL+1),2);
latvec_t1 = zeros((NL+1)*(NL+1),2);
for i1 = -NL:NL
    for i2 = -NL:NL
        for ias = 1:size(atomslc,1)
            atomssc_t1(tindx,:) = atomslc(ias,:) + [i1 i2];
            latvec_t1(tindx,:) = [i1 i2];
            tindx = tindx + 1;
        end
    end
end
latvec_t1 = int32(latvec_t1);

% new lattice vectors
% the second lattice vector is in the direction of the chain
% (NORMAL TO THE SURFACE)
avec2 = [1.0 0.0;
         0.0 1.0];
     
% get atoms in terms of the lattice coordinates of 
% the *new* lattice vectors
atomslc_t1 = inv(avec2)*atomssc_t1';
atomslc_t1 = atomslc_t1';

% *new* unit cell lattice atoms
inds = find((sum(abs(floor(atomslc_t1')))) < 1e-10);
atomslc2 = atomslc_t1(inds,:);

% *new* unite cell (in cartesian coords)
atomssc2 = avec2*atomslc2';
atomssc2 = atomssc2';

% determine NN atoms for "diag" -- (onsite?) 
%
% (0,0) lattice vector
%
tnnmat = [];
dist = 1.0;
for i = 1:size(atomssc2,1)
    for j = 1:size(atomssc2,1)
        v = atomssc2(j,:)-atomssc2(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat = [tnnmat; i, j];
        end
    end
end

% DEBUG
fprintf('onsite (0,0)\n');
for it = 1:size(tnnmat,1)
    i = tnnmat(it,1);
    j = tnnmat(it,2);
    dist = norm(atomssc2(j,:)-atomssc2(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2(j,:), atomssc2(i,:), dist);
end

% determine NN atoms for "off-diag" -- (offsite?) (DIRECTION NORMAL TO SURFACE)
% take the original atoms in the unit cell and add a lattice vector
% one unit cell away
%
% (0,1) lattice vector
%
atomssc2_t1 = atomssc2;
for i = 1:size(atomssc2_t1,1)
    atomssc2_t1(i,:) = atomssc2_t1(i,:) + avec2(:,2)';
end
tnnmat2 = [];
dist = 1.0;
for i = 1:size(atomssc2,1)
    for j = 1:size(atomssc2_t1,1)
        v = atomssc2_t1(j,:)-atomssc2(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat2 = [tnnmat2; i, j];
        end
    end
end

% DEBUG
fprintf('\n\n');
fprintf('offsite (0,1)\n');
for it = 1:size(tnnmat2,1)
    i = tnnmat2(it,1);
    j = tnnmat2(it,2);
    dist = norm(atomssc2_t1(j,:)-atomssc2(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2_t1(j,:), atomssc2(i,:), dist);
end

% determine NN atoms for "off-diag" -- (offsite?) (DIRECTION PARALLEL TO SURFACE)
% take the original atoms in the unit cell and add a lattice vector
% one unit cell away
%
% (1,0) lattice vector
%
atomssc2_t2 = atomssc2;
dist = 1.0;
for i = 1:size(atomssc2_t2,1)
    atomssc2_t2(i,:) = atomssc2_t2(i,:) + avec2(:,1)';
end
tnnmat3 = [];
dist = 1.0;
for i = 1:size(atomssc2,1)
    for j = 1:size(atomssc2_t2,1)
        v = atomssc2_t2(j,:)-atomssc2(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat3 = [tnnmat3; i, j];
            tnnmat3 = [tnnmat3; j, i];
        end
    end
end

% DEBUG
fprintf('\n\n');
fprintf('offsite (1,0)\n');
for it = 1:size(tnnmat3,1)
    i = tnnmat3(it,1);
    j = tnnmat3(it,2);
    dist = norm(atomssc2_t2(j,:)-atomssc2(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2_t2(j,:), atomssc2(i,:), dist);
end

% determine NN atoms for "off-diag" -- (offsite?) (DIRECTION ANGLE TO SURFACE)
% take the original atoms in the unit cell and add a lattice vector
% two unit cell away
%
% (1,1) lattice vector
%
atomssc2_t3 = atomssc2;
dist = 1.0;
for i = 1:size(atomssc2_t3,1)
    atomssc2_t3(i,:) = atomssc2_t3(i,:) + avec2(:,1)' + avec2(:,2)';
end
tnnmat4 = [];
dist = 1.0;
for i = 1:size(atomssc2_t1,1)
    for j = 1:size(atomssc2_t3,1)
        v = atomssc2_t3(j,:)-atomssc2_t1(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat4 = [tnnmat4; i, j];
        end
    end
end

% DEBUG
fprintf('\n\n');
fprintf('offsite (1,1)\n');
for it = 1:size(tnnmat4,1)
    i = tnnmat4(it,1);
    j = tnnmat4(it,2);
    dist = norm(atomssc2_t3(j,:)-atomssc2_t1(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2_t3(j,:), atomssc2_t1(i,:), dist);
end

% determine NN atoms for "off-diag" -- (offsite?) (DIRECTION PARALLEL TO SURFACE)
% take the original atoms in the unit cell and add a lattice vector
% one unit cell away
%
% (-1,0) lattice vector
%
atomssc2_t5 = atomssc2;
dist = 1.0;
for i = 1:size(atomssc2_t5,1)
    atomssc2_t5(i,:) = atomssc2_t5(i,:) - avec2(:,1)';
end
tnnmat5 = [];
dist = 1.0;
for i = 1:size(atomssc2_t1,1)
    for j = 1:size(atomssc2_t5,1)
        v = atomssc2_t5(j,:)-atomssc2_t1(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat5 = [tnnmat5; i, j];
            tnnmat5 = [tnnmat5; j, i];
        end
    end
end

% DEBUG
fprintf('\n\n');
fprintf('offsite (-1,0)\n');
for it = 1:size(tnnmat5,1)
    i = tnnmat5(it,1);
    j = tnnmat5(it,2);
    dist = norm(atomssc2_t5(j,:)-atomssc2_t1(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2_t5(j,:), atomssc2_t1(i,:), dist);
end

% determine NN atoms for "off-diag" -- (offsite?) (DIRECTION ANGLE TO SURFACE)
% take the original atoms in the unit cell and add a lattice vector
% two unit cell away
%
% (-1,1) lattice vector
%
atomssc2_t6 = atomssc2;
dist = 1.0;
for i = 1:size(atomssc2_t6,1)
    atomssc2_t6(i,:) = atomssc2_t6(i,:) + avec2(:,1)' + avec2(:,2)';
end
tnnmat6 = [];
dist = 1.0;
for i = 1:size(atomssc2,1)
    for j = 1:size(atomssc2_t6,1)
        v = atomssc2_t6(j,:)-atomssc2(i,:);
        if (abs(norm(v)-dist) < 1e-8)
            tnnmat6 = [tnnmat6; i, j];
        end
    end
end

% DEBUG
fprintf('\n\n');
fprintf('offsite (-1,1)\n');
for it = 1:size(tnnmat6,1)
    i = tnnmat6(it,1);
    j = tnnmat6(it,2);
    dist = norm(atomssc2_t6(j,:)-atomssc2(i,:));
    fprintf('it: %d    (%6.3f,%6.3f) - (%6.3f,%6.3f)  dist: %6.3f\n', ...
        it, atomssc2_t6(j,:), atomssc2(i,:), dist);
end

% build matrices
nsites = size(atomslc2,1);
M11_0 = zeros(nsites,nsites);
M11_1 = zeros(nsites,nsites);
M11_2 = zeros(nsites,nsites);
M12_0 = zeros(nsites,nsites);
M12_1 = zeros(nsites,nsites);
M12_2 = zeros(nsites,nsites);
M_offdiag = zeros(nsites,nsites);
for i = 1:size(tnnmat,1)
    M11_0(tnnmat(i,1),tnnmat(i,2)) = -thop;
end
for i = 1:size(tnnmat3,1)
    M11_1(tnnmat3(i,1),tnnmat3(i,2)) = -thop;
end
for i = 1:size(tnnmat5,1)
    M11_2(tnnmat5(i,1),tnnmat5(i,2)) = -thop;
end
for i = 1:size(tnnmat2,1)
    M12_0(tnnmat2(i,1),tnnmat2(i,2)) = -thop;
end
for i = 1:size(tnnmat4,1)
    M12_1(tnnmat4(i,1),tnnmat4(i,2)) = -thop;
end
for i = 1:size(tnnmat6,1)
    M12_2(tnnmat6(i,1),tnnmat6(i,2)) = -thop;
end

