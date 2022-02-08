% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel F, Diekmann O, Vermiglio R (2021). Numerical bifurcation analysis of renewal equations via pseudospectral approximation, 
% Journal of Computational and Applied Mathematics, 397:113611. https://doi.org/10.1016/j.cam.2021.113611

% MC_example22
% command line instructions for the MatCont continuation of the system defined
% in PS_example22.m
% (Gurney et al, Nature 1980)

clear;
clearvars -global cds
close all

%% Analytic Hopf curve

w = [1.7:0.01:3]'; % linspace(pi/2,pi-0.1,100)';

c1= w.*cos(w)./sin(w);
c2= -w./sin(w);

mmu = -c1;
bbeta = -c1.*exp(1+c2./c1);

figure
plot(mmu,bbeta./mmu); hold on
axis([0 5 0 50])
xlabel('mu'); ylabel('beta/mu')

%% Numerical continuation with Matcont

% Discretization parameters
M=20;

% Initial parameter values
tau_max=5; 
mu=4;
gamma=1;
aux=1;
par=[gamma,mu,aux,tau_max,M]';

% Approximated equilibrium corresponding to par
xeq=0;
yeq=[];

% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
ap2=2;
ap3=3;
TOL=1e-6;
TestTOL=1e-6;

%% Continuation process

MM=M; % dimension of the approximating ODE system
handles=feval(@PS_example22);
opt=contset;
global cds;

%% Equilibrium continuation from initial point [xeq;yeq]

display('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TestTOL);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxNumPoints',200);

state_eq=feval(handles{1},M,xeq,yeq,tau_max); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_example22,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=0;
while ((length(se)<3) && xe(end,end)>0.1 && jj<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
    jj=jj+1;
end

%% Detection of singular points
% xe,ve,se,he,fe = output of previous continuation
% par = current parameter vector
% ap1 = index of continuation parameter in vector par

% BP, branching point
for ii=1:length(se)
    if strcmp(se(ii).label,'BP')==1
        BP_index=se(ii).index;
        sBP=se(ii);
        break;
    end
end
par(ap1)=xe(end,BP_index);
BP=xe(1:MM,BP_index);

xeBP=xe; veBP=ve; seBP=se; heBP=he; feBP=fe;
parBP=par;

%% Equilibrium continuation from BP
% BP = vector of variables at BP
% parBP = parameter vector at BP
% sBP = information about BP from previous continuation
display('Starting equilibrium continuation from BP');

% set options
%opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'MaxStepsize',1);
opt=contset(opt,'Backward',0);

[x0,v0]=init_BP_EP(@PS_example22,BP,parBP,sBP,0.01);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=0;
while ( xe(end,end)<30 && xe(end,end)>0 && jj<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
    jj=jj+1;
end

% Plot
figure(1) 
cpl(xe,ve,se,[MM+1 1]);
xlabel('gamma');
title(['Bifurcation Example 2.2, mu=',num2str(mu),', M=',num2str(M)]);

%% Detection of singular points
% H, Hopf point

for ii=1:size(se)
    if (strcmp(se(ii).label,'H ')==1 && xe(end,se(ii).index)>0)
        H_index=se(ii).index;
        break;
    end
end
par(ap1)=xe(end,H_index);
H=xe(1:MM,H_index);

xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
parH=par;

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

% set options
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Singularities',0);

% forward branch
opt=contset(opt,'Backward',0);

[x0,v0]=init_H_H(@PS_example22,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)

jj=0;
while (xh(MM+2,1)<5 && xh(MM+2,end)<5 && xh(MM+2,1)>0.6 && xh(MM+2,end)>0.6 && jj<5)
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+2,end)
      jj=jj+1;
end

%% backward branch
opt=contset(opt,'Backward',1);

[x0,v0]=init_H_H(@PS_example22,H,parH,[ap1 ap2]);
[xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
% [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)

jj=0;
while (xhb(MM+2,1)<5 && xhb(MM+2,end)<5 && xhb(MM+2,1)>0.6 && xhb(MM+2,end)>0.6 && jj<5)
    [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+2,end)
    jj=jj+1;
end

%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

% set options
TOL=1e-6;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

opt=contset(opt,'MaxNumPoints',200);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepsize',10); %100
opt=contset(opt,'InitStepsize',1e-1);
opt=contset(opt,'Adapt',1);

ntst=40; % number of intervals
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_example22,H,parH,ap1,1e-2,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
jj=1;
while ((length(slc)<5) && xlc(end,end)< 90 && jj<10)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
    jj=jj+1;
end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

%% detection of PD bifurcation

SPD=1;
for S=1:size(slc)
    if strcmp(slc(S).label,'PD ')==1
        PD_index=slc(S).index;
        SPD=S;
        break;
    end
end

par_PD=xlc(end,PD_index);
parPD=parH;
parPD(ap1)=par_PD;

%% PD continuation in two parameters % ap1,ap2 = index of continuation
% parameters in the vector par display('Starting PD continuation');
% 
TOL=1e-6;
FunTOL=1e-6;

% set options
opt=contset(opt,'FunTolerance',FunTOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL); %opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Multipliers',0); 
opt=contset(opt,'Singularities',0);
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'MaxStepsize',100);

% forward branch
opt=contset(opt,'Backward',0);

[x0,v0]=init_PD_PD(@PS_example22,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
[xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); % 
l=size(xpd,1);

jj=0; 
while (xpd(l,1)<5 && xpd(l,end)<5 && xpd(l,1)>0.6 && xpd(l,end)>0.6 && jj<5) %
    [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); 
    xpd(l,end) %   
    xpd(end-1,end)./xpd(end,end)
    jj=jj+1;
end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

%% backward branch
opt=contset(opt,'Backward',1);
opt=contset(opt,'MaxStepsize',100);

[x0,v0]=init_PD_PD(@PS_example22,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
[xpd1,vpd1,spd1,hpd1,fpd1]=cont(@perioddoubling,x0,v0,opt); % 

jj=0; 
while (xpd1(l,1)<5 && xpd1(l,end)<5 && xpd1(l,1)>0.6 && xpd1(l,end)>0.6 && jj<5) %  
    [xpd1,vpd1,spd1,hpd1,fpd1]=cont(xpd1,vpd1,spd1,hpd1,fpd1,cds); xpd1(l,end) %       
    jj=jj+1;
end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

% PLOTS
%% Plot of stability regions

figure(2); clf
plot([0;5], [1;1], 'g'); hold on % existence
plot(xh(MM+2,:),xh(MM+1,:)./xh(MM+2,:),'b'); % Hopf
plot(xhb(MM+2,:),xhb(MM+1,:)./xhb(MM+2,:), 'b')
plot(xpd(end,:),xpd(end-1,:)./xpd(end,:),'r'); % Period doubling
plot(xpd1(end,:),xpd1(end-1,:)./xpd1(end,:),'r'); % Period doubling
axis([0 5 0 50]);
xlabel('mu'); ylabel('gamma/mu');
title(['Stability regions, M=',num2str(M),', tau=',num2str(tau_max)]);
% savefig(['M',num2str(M),'_PD_curve_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]); % ]); %

%% Plot of bifurcation diagram of first component
p = pi*(2*(0:M-1)'+1)/(2*M);
x=[1;sin(pi/2-p)]; % nodes with addition of 1 % either cos(p) or sin (pi/2-p) 
X=repmat(x,1,M+1);
dX=X-X';
c=[2^(M-1)/M*prod(dX(1,2:end)); ((-1).^(0:M-1)').*dX(2:end,1)./sin(p)];
D=(c*(1./c'))./(dX+(eye(M+1)));
D=D-diag(sum(D')); % differentiation matrix
% scaling
Nodes = 0.5*tau_max*(x-1);
DM = 2/tau_max*D(2:end,2:end);

[QuadWeights,QuadNodes]=cheb_quad(50,-tau_max,-1);

b0BP = zeros(size(xeBP,2),1);
for index_sol = 1:size(xeBP,2)
    
    mub = xeBP(end,index_sol);
    % Parameters and functions
    beta0 = gamma*exp(mub);
    
    derb = DM*xeBP(1:end-1,index_sol);
    der = polint(Nodes(2:end),derb,QuadNodes); % polint from Weideman-Reddy suite

    ARG = QuadWeights*(der.*exp(mub*QuadNodes));

    FM = beta0*ARG.*exp(-100*ARG); % scaling for ease of computation: b(t)=100*btilde
    b0BP(index_sol) = FM; % select component in node -tau
end

b0H = zeros(size(xeH,2),1);
for index_sol = 1:size(xeH,2)
    
    mub = xeH(end,index_sol);
    % Parameters and functions
    beta0 = gamma*exp(mub);
    
    derb = DM*xeH(1:end-1,index_sol);
    der = polint(Nodes(2:end),derb,QuadNodes); % polint from Weideman-Reddy suite

    ARG = QuadWeights*(der.*exp(mub*QuadNodes));

    FM = beta0*ARG.*exp(-100*ARG); % scaling for ease of computation: b(t)=100*btilde
    b0H(index_sol) = FM; % select component in node -tau
end

figure(10); clf; 
xlabel('$\mu$','Interpreter','latex');
title(['Bifurcation Example 2.2, mu=',num2str(mu),', M=',num2str(M)]);

plot(xeBP(end,:),b0BP,'b');
hold on
plot(xeH(end,:),b0H,'b');

for ii=2:length(seBP)-1
    index=seBP(ii).index;
    plot(xeBP(end,index),b0BP(index),'or');
end

for ii=2:length(seH)-1
    index=seH(ii).index;
    plot(xeH(end,index),b0H(index),'or');
end

Per_Solutions = zeros(ntst*ncol+1,size(xlc,2));
for ind_persol=1:size(xlc,2)
    for ind_mesh=1:ntst*ncol+1

        mub = xlc(end,ind_persol);
        % Parameters and functions
        beta0 = gamma*exp(mub);

        derb = DM*xlc((ind_mesh-1)*MM+1:(ind_mesh-1)*MM+MM,ind_persol);
        der_quad = polint(Nodes(2:end),derb,QuadNodes);
        
        ARG = QuadWeights*(der_quad.*exp(mub*QuadNodes));

        FM = beta0*ARG.*exp(-100*ARG); % scaling for ease of computation: b(t)=100*btilde
        b0_per = FM;
        Per_Solutions(ind_mesh,ind_persol) = b0_per;

    end
        
end
    
upperbound=max(Per_Solutions);
lowerbound=min(Per_Solutions);

plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'or',xlc(end,index),lowerbound(index),'or');
end
% savefig(['M',num2str(M),'_bif_diagram_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);


function [w,x]=cheb_quad(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% see Trefethen 2000

    p=pi*(0:N)'/N;
    x=((a-b)*cos(p)+b+a)/2;

    % Quadrature weights
    w=zeros(1,N+1);
    ii=2:N;
    v=ones(N-1,1);
    if mod(N,2)==0
        w(1)=1/(N^2-1);
        w(N+1)=w(1);
        for k=1:N/2-1
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
        v=v-cos(N*p(ii))/(N^2-1);
    else
        w(1)=1/N^2;
        w(N+1)=w(1);
        for k=1:(N-1)/2
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
    end
    w(ii)=2*v/N;
    w=w*abs(b-a)/2;

end

%%

function [w,x]=cheb_quad(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% see Trefethen 2000

    p=pi*(0:N)'/N;
    x=((a-b)*cos(p)+b+a)/2;

    % Quadrature weights
    w=zeros(1,N+1);
    ii=2:N;
    v=ones(N-1,1);
    if mod(N,2)==0
        w(1)=1/(N^2-1);
        w(N+1)=w(1);
        for k=1:N/2-1
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
        v=v-cos(N*p(ii))/(N^2-1);
    else
        w(1)=1/N^2;
        w(N+1)=w(1);
        for k=1:(N-1)/2
            v=v-2*cos(2*k*p(ii))/(4*k^2-1);
        end
    end
    w(ii)=2*v/N;
    w=w*abs(b-a)/2;

end
