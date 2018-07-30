function [zout,tout,fout,dn,rout,dcon,tcon]=solvestagnantrift(vn,Di,etat,k0,K0,eta0,H0,Tb,Ts,Fb,kfun,eta,hfun,Hs,DT,Ra,Th,T0,L0,F0) %% 

%% Obtain convective core;
[bcv,DTv,Ttv,Tiv,Hiv,Ftv,bbv,btv,biv,etai,etav]=...
    convectivecore(Di,etat,k0,K0,eta0,H0,Tb,Fb,kfun,eta,hfun,Hs,DT,Ra,Th);
Rcv=Tiv*bcv+(Ttv-Tiv)*btv/2+(Tb-Tiv)*bbv/2-Tb*bcv;

%% Add conductive layer
% assume no upwelling
% Scales;
kv=k0*kfun(Ttv/Tb);
Kv=K0*kfun(Ttv/Tb);
Hv=H0*hfun(Ttv/Tb);

%% Scales
Tv=Ttv-Ts;
Lv=sqrt(kv*Tv/Hv);
Vv=Kv/Lv;
Fv=Hv*Lv;

fb=Ftv/Fv;
dz=1e-4;
nz=100; 
tb=Ttv/Tv;           % bottom temperature (ND)
% rlv0=rhoi*Lh*V0;    % freezing

[znv,tnv,fnv,kout,hout,dnv,rnv]=...
        ConductionSolveRift(dz,tb,fb,vn,kfun,hfun,nz,Di/L0);
Dnv=dnv*Lv;
Znv=(dnv-znv)*Lv;
Tnv=tnv*Tv;
Fnv=fnv*Fv;
Rnv=rnv*Tv*Lv;


%% Merge
dn=(Dnv+bcv)/L0; %total thickness
zout=([Dnv+[bcv,btv+biv,btv],Znv])/L0; %depth vector
tout=[Tb,Tiv,Tiv,Tnv]/T0; %ND temperature vector
fout=[Fb,Fb,Ftv,Fnv]/F0; %ND flux
rout=(Rcv+Rnv)/T0/L0; % ND integral temperature
dcon=Dnv/L0;
tcon=Ttv/T0;
% kout=k0*kfun(tout/Tb);
% hout=H0*hfun(tout/Tb);
