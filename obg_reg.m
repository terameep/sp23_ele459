function [L,delta1,delta2]=obg_reg(A,B,C,K,poles,T)
%OBG_REG  Observer gain calculation for observer-based regulator
%[L,delta1,delta2] =obg_reg(A,B,C,K,poles,T)
%
%INPUTS
%  A,B,C        State-space plant model (continuous- or discrete-time)
%  K            State-feedback gain matrix
%  poles        Vector of desired observer poles
%  T            Sampling interval (Use T=0 for continuous-time system)
%
%OUTPUT
%  L            Observer gain matrix. Note: eig(A-L*C)=poles
%  delta1       Unstructured input-multiplicative infinity-norm 
%               perturbation bound 
%  delta2       Unstructured input-multiplicative feedback 
%               infinity-norm perturbation bound
%
%
% R.J. Vaccaro 3/2016
%
[m,n]=size(C);
cvo=~[imag(poles)==0];
ind=find(cvo);
if length(ind>1)
    ind=ind(1:2:end);
    cvo(ind)=[];
    poles(ind)=[];
end
[W1,W2]=obasis2(A',C',cvo,poles);
ncp=2*sum(cvo); % number of complex poles
nrp=n-ncp; % number of real poles
ncvo=length(cvo);
nt=n*(m-1); % number of theta parameters needed
options=optimset('display','off','MaxFunEvals',5000);
theta=pi/4*ones(nt,1);
%%%%%%%%%%%%%%%%%%%%

tfun=[.1 .01];
tfun=[tfun ones(1,6)*1.e-5];
tx=[0.01 0.001];
tx=[tx ones(1,6)*1.e-7];
tfun=[.01 .01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
%tfun=[tfun ones(1,1)*1.e-5];
tx=[0.01 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001];
fprintf('\ndelta1   delta2 \n')
for LOOP=1:4
    options=optimset('display','off','TolX',tx(LOOP),'TolFun',tfun(LOOP),...
        'MaxFunEvals',5000);
[theta,f]=fminsearch(@cost,theta,options,A,B,C,K,cvo,W1,W2,T);
L=formL(W1,W2,cvo,theta,n,m)';
[delta1,delta2]=rb_regob(A,B,C,K,L,T);
fprintf('%1.4f   %1.4f\n',delta1,delta2)
end
L=formL(W1,W2,cvo,theta,n,m)';

%
function [L,cond_flg]=formL(V1,V2,cv,theta,n,p)
cond_flg=0;
ncp=sum(cv); % number of complex poles/2
nrp=n-2*ncp; % number of real poles
M=[];
if nrp>0
    for k=1:nrp
        theta1=theta((k-1)*(p-1)+1:k*(p-1));
        alpha1=hyperspherical(theta1);
        M=[M V1(:,(k-1)*p+1:k*p)*alpha1];
    end
end
if ncp>0
    offset1=nrp*(p-1);
    offset2=offset1+ncp*(p-1);
    for k=1:ncp
        theta1=theta(offset1+(k-1)*(p-1)+1:offset1+k*(p-1));
        alpha1=hyperspherical(theta1);
        theta2=theta(offset2+(k-1)*(p-1)+1:offset2+k*(p-1));
        alpha2=alpha1.*[1;exp(j*theta2)];
        X=V2(:,(k-1)*p+1:k*p)*alpha2;
        M=[M real(X) imag(X)];
    end
end
if cond(M(1:n,:))>1.e9
    L=[];
    cond_flg=1;
else
    L=M(n+1:end,:)/M(1:n,:);
end

function h=hyperspherical(theta)
%theta is a vector containing n-1 angles (radians)
%h is a vector containing n Cartesian coordinates of
%a point on a unit hypersphere
p=length(theta);
h=[cos(theta);1];
for k=1:p
    h=h.*[ones(k,1);sin(theta(k))*ones(p-k+1,1)];
end

%
function f=cost(theta,A,B,C,K,cvo,V1,V2,T)
[m,n]=size(C);
[L,cond_flg]=formL(V1,V2,cvo,theta,n,m);
L=L';
if ~cond_flg
[n,p]=size(B);
D1=zeros(p,p);
    sys1=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],D1,T);
    f1=norm(sys1,inf);
    sys2=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],eye(p),T);
    f2=norm(sys2,inf);
    f=f1+100*f2;
else
    f=1.e10;
end
%
function [V1,V2]=obasis2(A,B,cv,poles)
[n,p]=size(B);
ncv=length(cv);
scv=sum(cv);
V1=[];
V2=[];
if n-2*scv>0
    V1=zeros(n+p,(n-2*scv)*p); % basis vectors for real-valued poles
end
if scv>0
    V2=zeros(n+p,scv*p); % basis vectors for complex-valued poles
end
c1=0;
c2=0;
for k=1:length(cv)
    X=null([(poles(k)*eye(n)-A) B]);
    if ~cv(k)
        c2=c2+1;
        V1(:,(c2-1)*p+1:c2*p)=X;
    else
        c1=c1+1;
        V2(:,p*(c1-1)+1:p*c1)=X;
    end
end
function [delta1,delta2]=rb_regob(A,B,C,K,L,T)
%RB_REGOB       Computes input-multiplicative stability robustness bound
%               for observer-based regulator
% [delta1,delta2]=rb_regob(A,B,C,K,L,T)
%INPUTS
%  A,B,C        State-space plant model (continuous- or discrete-time)
%  K            State feedback gain matrix
%  L            Observer gain matrix
%  T            Sampling interval (Use T=0 for continuous-time system)
%
%OUTPUT
%  delta1       Unstructured input-multiplicative infinity-norm 
%               perturbation bound 
%  delta2       Unstructured input-multiplicative feedback 
%               infinity-norm perturbation bound
%
% R.J. Vaccaro 3/2016
[n,p]=size(B);
D1=zeros(p,p);
    sys1=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],D1,T);
    delta1=1/norm(sys1,inf);
    sys2=ss([A -B*K;L*C A-L*C-B*K],[B; 0*B],[0*K -K],eye(p),T);
    delta2=1/norm(sys2,inf);

