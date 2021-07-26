function out = PS_example22
% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel F, Diekmann O, Vermiglio R (2021). Numerical bifurcation analysis of renewal equations via pseudospectral approximation, 
% Journal of Computational and Applied Mathematics, 397:113611. https://doi.org/10.1016/j.cam.2021.113611

%% PS_example22.m
% MatCont system definition file of the PseudoSpectral Discretization of
% the nonlinear RE defined in the paper Gurney et al, Nature, 1980,
% x(t) = beta*h(int_1^amax x(t-a)*exp(x(t-a))da)
% with h(x)= x*exp(-x)
% for the integrated state B=int_0^t x(s)ds
% using Chebyshev zeros (plus 0)
% The code uses the function polint.m, available from the Differentiation
% Matrix Suite by Weideman, Reddy, 2000

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= @userf; % user function to select specific parameter values 
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,gamma,mu,aux,tau,M) 

%% discretization of the unitary interval [-1,0]
    % construction of nodes and differentiation matrix
    p = pi*(2*(0:M-1)'+1)/(2*M);
    x=[1;sin(pi/2-p)]; % nodes with addition of 1 % either cos(p) or sin (pi/2-p) 
    X=repmat(x,1,M+1);
    dX=X-X';
    c=[2^(M-1)/M*prod(dX(1,2:end)); ((-1).^(0:M-1)').*dX(2:end,1)./sin(p)];
    D=(c*(1./c'))./(dX+(eye(M+1)));
    D=D-diag(sum(D')); % differentiation matrix

    % scaling
    Nodes = 0.5*tau*(x-1);
    DD = 2/tau*D;
    DM = DD(2:end,2:end);

    %% SYSTEM DEFINITION *** specific to the equation ***

    % Parameters and functions
    abar=1;
    beta0 = gamma*exp(mu);

    % For quadrature formulas 
    [QuadWeights,QuadNodes]=cheb_quad(50,-tau,-abar);

    der_state = DM*state;
    der = polint(Nodes(2:end),der_state,QuadNodes); % polint from Weideman-Reddy suite

    ARG = QuadWeights*(der.*exp(mu*QuadNodes));
    % ARG = integral(@(theta) polint(tau*UnitNodes(2:end),der_state,theta)'.*exp(mu*theta),-tau,-abar);

    % FM = beta0*ARG.*exp(-ARG);
    FM = beta0*ARG.*exp(-100*ARG); % scaling for ease of computation: b(t)=100*btilde

    %% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

    dydt= der_state - FM;

end
 
 
% --------------------------------------------------------------------------
function Weq=init(M,xeq,yeq,tau)
% INPUT: M is the discretization parameter
%        xeq,yeq are column vectors with the equilibrium states of the RE and
%        DDE respectively
% OUTPUT Weq is the initial vector for init_EP_EP
    
    [~,UnitNodes]=cheb_quad(M,0,-1);
    Nodes=tau*UnitNodes;
    Weq=[kron(Nodes(2:end),xeq); kron(ones(M+1,1),yeq)];
    
end


function out=userf(time,state,gamma,mu,aux,tau,M) 
% Userfunction to select specific value of parameter
    
    out = (gamma-60).*(gamma-70);
    
end

% ------------
%%% AUXILIARY FUNCTIONS

function [w,x]=cheb_quad(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
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
