%% Fast Hopf

clear all
close all

NPARCELLS=360;
load('avg_SC.mat')
C = SC_avg_weights;
C = C(1:360,1:360);
SC_avg_dists(C < prctile(C(:),50)) = 0;
C(C < prctile(C(:),50)) = 0; % threshold SC
C = C./ 4e5; % normalize SC
C = C.^(1/2);
C(C>0) = C(C>0) - min(C(C>0));

% optimal working point
a = -0.008;
G = 0.22025;

% Parameters HOPF
Tmax=1200;
TR=1;
f_diff=0.055*ones(1,NPARCELLS);
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.02;
dsig = sqrt(dt)*sig;


wC = G*squeeze(C);
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

%% Hopf Simulation

xs=zeros(Tmax,NPARCELLS);
z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
nn=0;
% discard first 2000 time steps
for t=0:dt:2000
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
end
% actual modeling
for t=0:dt:((Tmax-1)*TR)
    suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
    zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
    z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    if abs(mod(t,TR))<0.01
        nn=nn+1;
        xs(nn,:)=z(:,1)';
    end
end
