close all;
clear all;

% amplitude
A = 12.5;
% size of the box
L = 5.0;
% number of nstates
nstates = 2;
% max G vector length
maxglen = 4.0;

% compute full up 1-d cosine pw matrix
[C, kmesh, gvecs, e] = diag_cos_pw(A, L, maxglen, 1, nstates, 0, 0);
ngvecs = size(gvecs,1);

% get gmax and gmin
gmax = round(max(gvecs(:,1))*L/2/pi);
gmin = round(min(gvecs(:,1))*L/2/pi);
npts = gmax - gmin + 1;

% check number of states
if (nstates > npts)
    % would like npts to be odd [... -2 -1 0 1 2 ...]
    if (mod(nstates,2)==0)
        npts = nstates + 1;
    end
end

% print output for debug
gmax2 = round(max(gvecs*L/2/pi));
gmin2 = round(min(gvecs*L/2/pi));
fprintf('npts: %5d\n\n',npts);
fprintf('gmin: %7d%7d%7d\n',gmin2(1),gmin2(2),gmin2(3));
fprintf('gmax: %7d%7d%7d\n',gmax2(1),gmax2(2),gmax2(3));
fprintf('\n\n');

% generate 1d g-vectors
gvecs1d = [-(npts-1)/2:(npts-1)/2]'*2*pi/L;

% test code
% ww = [gvecs1d abs(gvecs1d)];
% ww = sortrows(ww,2);
% gvecs1d = ww(:,1);

% for z-component
% generate potential for -A(cos(2*pi*n*x/L) + 1)
% Vconst = -A*diag(ones(1,npts));
% Vcos = -0.5*(diag(ones(npts-1,1),1) + diag(ones(npts-1,1),-1));
% V = Vconst + Vcos;
% % generate T
% T = 0.5*diag(gvecs(:,3).^2);
% H = T + V;
%[cz,evHz] = eig(H);

% for x and y components
cx = eye(npts);
evHx = 0.5*(gvecs1d.^2);
cy = eye(npts);
evHy = 0.5*(gvecs1d.^2);
cz = eye(npts);
evHz = 0.5*(gvecs1d.^2);

% now we need to put the x, y and z components back together intelligently
CC = zeros(size(C));
igvecs = round(gvecs*L/2/pi);
igvecs1d = round(gvecs1d*L/2/pi);
igvecs1d2 = igvecs1d+min(igvecs1d)+1;

%for ig = 1:ngvecs
tmp_states = zeros(3000,3);
itst = 0;
for ig = 1:30
    % get current g-vector
    igvec = igvecs(ig,:);
    g1 = igvecs(ig,1); g2 = igvecs(ig,2); g3 = igvecs(ig,3);
    % what indicies into igvecs1d correspond to current g-vector?
    [t1, ig1] = min(abs(igvecs1d-g1));
    [t1, ig2] = min(abs(igvecs1d-g2));
    [t1, ig3] = min(abs(igvecs1d-g3));
    if ((abs(g1-igvecs1d(ig1)) > 1e-8) || ...
        (abs(g2-igvecs1d(ig2)) > 1e-8) || ...
        (abs(g3-igvecs1d(ig3)) > 1e-8))
        error('error: was not able to find g-vector\n');
    end
    for ist1 = 1:npts
        for ist2 = 1:npts
            for ist3 = 1:npts
                t2 = cx(ig1,ist1)*cy(ig2,ist2)*cz(ig3,ist3);
                if (abs(t2) > 1e-5)                    
                    t3 = evHx(ig1) + evHy(ig2) + evHz(ig3);
                    itst = itst + 1;
                    tmp_states(itst,:) = [ig t2 t3];
                    fprintf('%4d%4d%4d%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\n',...
                        g1,g2,g3,cx(ig1,ist1),...
                        cy(ig2,ist2),cz(ig3,ist3),evHx(ig1),...
                        evHy(ig2),evHz(ig3),t3);
                end
            end
        end
    end
    %CC(1,ig,ist) = cx(ig1,ist)*cy(ig2,ist)*cz(ig3,ist);
    end

% % loop over all states and gvectors and build CC
% for ist = 1:2
%     fprintf('\nist: %4d\n',ist);
%     %for ig = 1:ngvecs
%     for ig = 1:2
%         % get current g-vector
%         igvec = igvecs(ig,:);
%         g1 = igvecs(ig,1); g2 = igvecs(ig,2); g3 = igvecs(ig,3);
%         % what indicies into igvecs1d correspond to current g-vector?
%         [t1, ig1] = min(abs(igvecs1d-g1));
%         [t1, ig2] = min(abs(igvecs1d-g2));
%         [t1, ig3] = min(abs(igvecs1d-g3));
%         if ((abs(g1-igvecs1d(ig1)) > 1e-8) || ...
%             (abs(g2-igvecs1d(ig2)) > 1e-8) || ...
%             (abs(g3-igvecs1d(ig3)) > 1e-8))
%             error('error: was not able to find g-vector\n');
%         end
% %          CC(1,ig,ist) = cxy(ig1,ist)*cxy(ig2,ist)*cz(ig3,ist);
% %         fprintf('%15.8f%15.8f%15.8f\n',cxy(ig1,ist),cxy(ig2,ist),cz(ig3,ist));
%         CC(1,ig,ist) = cx(ig1,ist)*cy(ig2,ist)*cz(ig3,ist);
%         fprintf('%4d%4d%4d%15.8f%15.8f%15.8f\n',g1,g2,g3,cx(ig1,ist),...
%             cy(ig2,ist),cz(ig3,ist));
%     end
% end

% % build 3d states
% for ig1 = 1:npts
%     for ig2 = 1:npts
%         for ig3 = 1:npts
%             gvec = [gvecs1d(ig1) gvecs1d(ig2) gvecs1d(ig3)];
%             for ist1 = 1:npts
%                 for ist = 1:npts
%                     for ist3 = 1:npts
%                         t1 = cx(ig1,ist1)*cy(ig2,ist2)*cz(ig3,ist3);
%                         if (abs(t1) > 1e-6)
%                             
%                     end
%                 end
%             end
%             t1 = cx(
%         end
%     end
% end

