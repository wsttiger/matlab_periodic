
clear all;
close all;

% amplitude
A = 22.5;
% size of the box
L = 5.0;
% number of states
%nstates = 25;
nstates = 16;

% number of g-vectors
imax = 15;
%gmaxlen = 10.0;
gmaxlen = 10.0;

% generate g-vectors
fprintf('Generating g-vectors with length < %15.8f\n', gmaxlen);
gvecs = zeros(3000,3);
glens = zeros(3000,1);
ngvecs = 0;
t1 = 2*pi/L;
for i1 = -imax:imax
    for i2 = -imax:imax
        for i3 = -imax:imax
            glen = norm([i1 i2 i3]*t1);
            if (glen < gmaxlen)                
                ngvecs = ngvecs+1;
                gvecs(ngvecs,:) = [i1 i2 i3]*t1;
                glens(ngvecs) = glen;
            end
        end
    end
end

% sort g-vectors
gt1 = [gvecs(1:ngvecs,1:3), glens(1:ngvecs,1)];
gt1 = sortrows(gt1,4);
gvecs = gt1(1:ngvecs,1:3);
glens = gt1(1:ngvecs,4);

% generate potential for -A*(cos(2*pi*n*x/L)*cos(2*pi*n*y/L)*cos(2*pi*n*z/L) + 1)
fprintf('Generating potential for A = %15.8f\n', A);
fprintf('with boxsize = %15.8f\n\n', L);
tic
Vconst = -A*diag(ones(1,ngvecs));
Vcos = zeros(ngvecs,ngvecs);
% for the cosine potential we a series of 8 delta functions to search
% 8 = 2^3
% and what we're really searching for of the condition of V(G-G') where we
% are searching for G'
qvecs = zeros(8,3);
qvecs(1,:) = [ 1  1  1];
qvecs(2,:) = [ 1  1 -1];
qvecs(3,:) = [ 1 -1  1];
qvecs(4,:) = [ 1 -1 -1];
qvecs(5,:) = [-1  1  1];
qvecs(6,:) = [-1  1 -1];
qvecs(7,:) = [-1 -1  1];
qvecs(8,:) = [-1 -1 -1];
qvecs = qvecs*t1;
% loop over G
ghits = zeros(20000,3);
ghitcnt = 0;
for ig = 1:ngvecs
    % loop over G'
    % loop over qvecs
    for iqv = 1:8
        gt2 = gvecs(ig,:) - qvecs(iqv,:);
        for ig2 = 1:ngvecs
            if (norm(gvecs(ig2,:)-gt2) < 1e-8)
                Vcos(ig,ig2) = -0.125*A;
                ghitcnt = ghitcnt + 1;
                ghits(ghitcnt,:) = [ig ig2 iqv];
            end
        end
    end
end
vtime = toc;
fprintf('Time to make potential: %10.4f\n\n', vtime);

% debug info
ghitr1 = gvecs(ghits(1:ghitcnt,1),:) ./ t1;
ghitr2 = gvecs(ghits(1:ghitcnt,2),:) ./ t1;
ghitdf = ghitr1-ghitr2;

% k-points
% kpoints = [0.0, 0.0, 0.0;
%            0.0, 0.0, 0.5;
%            0.0, 0.5, 0.0;
%            0.0, 0.5, 0.5;
%            0.5, 0.0, 0.0;
%            0.5, 0.0, 0.5;
%            0.5, 0.5, 0.0;
%            0.5, 0.5, 0.5];
% kpoints = [0.0, 0.0, 0.0;
%            0.5, 0.5, 0.5];
kpoints = [0.5, 0.5, 0.5];
kpoints = kpoints*t1;
nkpts = size(kpoints,1);

% sparse
Vconst_S = sparse(Vconst);
Vcos_S = sparse(Vcos);
clear Vconst;
clear Vcos;

% file info
fname = 'basis';
fid = fopen(fname,'w');
fprintf(fid,'%d\n',ngvecs);
for i = 1:ngvecs
    fprintf(fid,'%15.8f    %15.8f    %15.8f\n', gvecs(i,1), gvecs(i,2), gvecs(i,3)); 
end
fprintf(fid,'%d\n',nkpts);
fprintf(fid,'%d\n',nstates);

Es = zeros(nstates,nkpts);
for kindx = 1:nkpts
% generate the kinetic energy matrix
    kvecs = [kpoints(kindx,1)*ones(ngvecs,1) kpoints(kindx,2)*ones(ngvecs,1) kpoints(kindx,3)*ones(ngvecs,1)];
    gkvecs = gvecs + kvecs;
    gkvecsnorms = zeros(ngvecs,1);
    for i = 1:ngvecs
    gkvecsnorms(i) = norm(gkvecs(i,:));
    end
    T = 0.5*diag(gkvecsnorms.^2);
    T_S = sparse(T);
    clear T;

    % full hamiltonian
    H_S = (T_S + Vconst_S + Vcos_S);

    [cs,es] = eigs(H_S,nstates,'sa');
    for ig = 1:ngvecs
        for is = 1:nstates
            fprintf(fid,'%15.8f    ', cs(ig,is));
        end
    fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fclose(fid);