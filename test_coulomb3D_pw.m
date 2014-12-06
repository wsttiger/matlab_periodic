
clear all;
close all;

% amplitude
A = 2.5;
% size of the box
L = 5.0;
% number of states
nstates = 5;

% number of g-vectors
imax = 15;
gmaxlen = 10.0;


% generate g-vectors
fprintf('Generating g-vectors with length < %15.8f\n', gmaxlen);
gvecs = zeros(1000,3);
glens = zeros(1000,1);
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
Vconst = -A*L*L*L*diag(ones(1,ngvecs));
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
                Vcos(ig,ig2) = -0.125*A*L*L*L;
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

% generate the kinetic energy matrix
fprintf('Kinetic energy matrix\n');
tic
T = 0.5*L*L*L*diag(glens.^2);
ttime = toc;
fprintf('Time to make KE matrix: %10.4f\n\n', ttime);

% full hamiltonian
H = (T + Vconst + Vcos)/L/L/L;
V = Vconst + Vcos;

Hs = sparse(H);
[Cs,Es] = eigs(Hs,nstates,'sa');
%Es = Es./L/L/L;
dEs = diag(Es);
