clear all;
close all;

%% Lattice constant in atomic units
a = 1.0; 
%% So we're assuming a gaussian potential of the form V = Bexp(-px^2)
%% where B will be a negative number
%% Amplitude of the potential (must be negative)
B = -120.0;
%% Parameter in the gaussian potential
p = 2.0;
%% Normalize to B
C = B * sqrt(p/pi);
%% Create G vectors
gtop = 20;
gsize = 2*gtop + 1;
G = 2*pi/a * [-gtop:gtop]';
% gsize = 100;
% G = 2*pi/a * [0:gsize-1];
%% Create wavevectors and energies
nkpts = 200;
kpoints = linspace(-pi/a, pi/a, nkpts);
%% ground and excited states of KP
en = zeros(size(kpoints));
ex = zeros(size(kpoints));
%% free electron gas energies
enfe = zeros(size(kpoints));
%% ground state wavefunctions of KP
phis = zeros(gsize,nkpts);
%% Loop over wavevectors
it = 1;
for k = kpoints
    %% Create matrix to be diagonalized
    A = zeros(gsize,gsize);
    FE = zeros(gsize,gsize);
     for gi = 1:gsize
        for gj = 1:gi-1
            %% Calculate fourier transform of potential Vf(Gi - Gj)
            Vgij = (C/sqrt(2*p))*exp(-((G(gi) - G(gj)).^2)/(4*p));
            A(gi,gj) = Vgij;
            A(gj,gi) = Vgij;
        end
       A(gi,gi) = C/sqrt(2*p) + 0.5*(k + G(gi)).^2;
       FE(gi,gi) = 0.5*(k + G(gi)).^2;
     end

    %% Do diagonalization
    [tmp, energies] = eig(A);
    %% ground state energies and wavefunctions
    en(it) = energies(1,1);
    phis(:,it) = tmp(:,1);
    %% 1st excited states
    ex(it) = energies(2,2);
    %% free electron gas
    energies = eig(FE);
    enfe(it) = min(energies);
    it = it + 1;
end

%% Plot potential
nxpnts = nkpts;
xpoints = linspace(-a/2,a/2, nxpnts);
V = C*exp(-p*xpoints.^2);
figure()
plot(xpoints, V, 'r');
xlabel('x');
ylabel('V');
title('potential in real space');
%% Plot dispersion
figure()
hold on
plot(kpoints, en);
plot(kpoints, ex);
hold off
xlabel('k');
ylabel('E');
title('dispersion of Kronig-Penney Model');
figure()
plot(kpoints, enfe, 'g');
xlabel('k');
ylabel('E');
title('dispersion of free electron gas');
%% Debug a little bit
delx = xpoints(2) - xpoints(1);
normV = sum(V)*delx;



