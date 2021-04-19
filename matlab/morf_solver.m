function [x,z_b_sym, dt]=morf_solver(B, delta_z_init, t_end)
%% 1D morphology computation (Setting up the model equations)

[x_sol, h]=B.solve;
x=linspace(min(x_sol), max(x_sol),numel(x_sol));
dx=median(diff(x));
z_bed_level=B.bed_level;
[x_bed_level, idx]=sort([B.x0, B.x_end]);
z_bed_level=z_bed_level(idx);
z_b_init=interp1(x_bed_level, z_bed_level, x) + interp1(x_sol,delta_z_init,x);
C=B.Chez;
S0=B.So;

poros=0.4;%B.porosity;
a_sed=2e-4;
b_sed=5;
u=B.Q./B.b./h;
qs=a_sed.*u.^b_sed;
c_b=b_sed./(1-poros)*qs./h;

dt=0.9*dx/max(c_b);
t=0:dt:t_end;
nt=numel(t);

dqsdu=b_sed*a_sed*u.^(b_sed-1);
K=1/(1-poros)*dqsdu*C^2.*h./3./u;



nx=numel(x);
mdiags=zeros(nx*nt,8); % matrix containing bands of sparse matrix
bnd=zeros(nx*nt,1);

% Initial condition
mdiags(1:nx,6)=1; bnd(1:nx)=z_b_init;

% Boundary conditions
% idx=1+(1:nt-1)*nx; mdiags(idx,6)=1; bnd(idx)=z_b(1);  % upstream
% idx=  (2:nt  )*nx; mdiags(idx,6)=1; bnd(idx)=z_b(end);% downstream

%%%%%%%%%%%%%%%%     Downstream boundary (advection equation)
% temporal derivatives
idx=bsxfun(@plus,nx,(1:nt-1)*nx); idx=idx(:); 
mdiags(idx, 6)=mdiags(idx, 6)+1/dt;
mdiags(idx-nx, 2)=mdiags(idx-nx, 2)-1/dt;

% advective term (1st order upwind)
idx=nx*ones(nt-1,1);
A=c_b(idx)'/dx;
rhs=-c_b(idx)*S0;
idx=bsxfun(@plus,nx,(1:nt-1)*nx); idx=idx(:); 
mdiags(   -1+idx, 5)=mdiags(   -1+idx,5)-A; %i-1, n
mdiags(      idx, 6)=mdiags(     +idx,6)+A; %i  , n
bnd(idx)=rhs;

%%%%%%%%%%%%%%%%      Temporal derivative (Euler forward)
idx=bsxfun(@plus,(1:nx-1)',(1:nt-1)*nx); idx=idx(:); 
mdiags(idx, 6)=mdiags(idx, 6)+1/dt; % i, n
mdiags(idx-nx, 2)=mdiags(idx-nx, 2)-1/dt; % i, n-1

%%%%%%%%%%%%%%%       Diffusive term (Crank Nicolson)
idx=repmat(2:nx-1,1,nt-1); 
D=-K(idx)'/2/dx^2;

idx=reshape(bsxfun(@plus,(2:nx-1)'  ,(1:nt-1)*nx),[],1); 
mdiags(    1+idx, 7)=mdiags(    1+idx, 7)+  D; % i+i, n
mdiags(      idx, 6)=mdiags(      idx, 6)-2*D; % i  , n
mdiags(   -1+idx, 5)=mdiags(   -1+idx, 5)+  D; % i-1, n
mdiags(-nx+1+idx, 3)=mdiags(-nx+1+idx, 3)+  D; % i+1, n-1
mdiags(-nx  +idx, 2)=mdiags(-nx  +idx, 2)-2*D; % i  , n-1
mdiags(-nx-1+idx, 1)=mdiags(-nx+1+idx, 1)+  D; % i-1, n-1

           
%%%%%%%%%%%%%%%%      Advective term (Upwind)

% 1st order for downstream side - 1 
idx=repmat(nx-1,1,nt-1);
A=-K(idx)'./c_b(idx)'/dx/dt;
idx=reshape(nx-1+(1:nt-1)*nx,[],1); 
mdiags(      idx, 6)=mdiags(      idx,6)-A; %i, n
mdiags(    1+idx, 7)=mdiags(    1+idx,7)+A; %i+1, n
mdiags(-nx  +idx, 2)=mdiags(-nx  +idx,2)+A; %i, n-1
mdiags(-nx+1+idx, 3)=mdiags(-nx+1+idx,3)-A; %i+1, n-1

% 2nd order for upstream side
idx=ones(1,nt-1);
A=-K(idx)'./c_b(idx)'/2/dx/dt;
idx=reshape(1+(1:nt-1)*nx,[],1);
mdiags(      idx, 6)=mdiags(      idx,6)-3*A; %i, n
mdiags(    1+idx, 7)=mdiags(    1+idx,7)+4*A; %i+1, n
mdiags(    2+idx, 8)=mdiags(    2+idx,8)-  A; %i+2, n
mdiags(-nx  +idx, 2)=mdiags(-nx  +idx,2)+3*A; %i, n-1
mdiags(-nx+1+idx, 3)=mdiags(-nx+1+idx,3)-4*A; %i+1, n-1
mdiags(-nx+2+idx, 4)=mdiags(-nx+2+idx,4)+  A; %i+2, n-1

% 3rd order for the rest
idx=repmat(2:nx-2,1,nt-1);
A=-K(idx)'./c_b(idx)'/6/dx/dt;
idx=reshape(bsxfun(@plus,(2:nx-2)',(1:nt-1)*nx),[],1);
mdiags(   -1+idx, 5)=mdiags(   -1+idx,5)-2*A; %i-1, n
mdiags(     +idx, 6)=mdiags(     +idx,6)-3*A; %i  , n
mdiags(   +1+idx, 7)=mdiags(   +1+idx,7)+6*A; %i+1, n
mdiags(   +2+idx, 8)=mdiags(   +2+idx,8)-  A; %i+2, n
mdiags(-nx-1+idx, 1)=mdiags(-nx-1+idx,1)+2*A; %i-1, n-1
mdiags(-nx  +idx, 2)=mdiags(-nx  +idx,2)+3*A; %i  , n-1
mdiags(-nx+1+idx, 3)=mdiags(-nx+1+idx,3)-6*A; %i+1, n-1
mdiags(-nx+2+idx, 4)=mdiags(-nx+2+idx,4)+  A; %i+2, n-1

%% 1D morphological model (Solving the equations)
M=spdiags(mdiags,[-nx-1 -nx -nx+1 -nx+2 -1 0 1 2],nx*nt,nx*nt);
z_b_sym=reshape(M\bnd,numel(x),nt);
