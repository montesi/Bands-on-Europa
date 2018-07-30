function [bcv,DTv,Ttv,Tiv,Hiv,Ftv,bbv,btv,biv,etai,etav]=convectivecore(Di,etat,k0,K0,eta0,H0,Tb,Fb,kfun,eta,hfun,Hs,DT,Ra,Th)

% Di:   maximum thickness
%% Demonstration

nb=100;
ball=linspace(0,Di,nb);
DTall=NaN(size(ball));
RaT=DTall;Thall=DTall;Ttop=DTall;Tint=DTall;
Tbot=Tb*ones(size(ball));

for ib=1:nb;
    bn=ball(ib);
    if bn~=0;
        [DTall(ib),RaT(ib),Thall(ib),Ttop(ib),Tint(ib)]=...
            TiIterate(bn,k0,K0,eta0,H0,Tb,Fb,kfun,eta,hfun,DT,Ra,Th);
    end
    %disp(sprintf('At b=%g km, DT=%g K',bn/1000,DTall(ib)));
end
Hall=Hs(ball,H0*hfun(Tint/Tb),k0*kfun(Tint/Tb),DTall);
% Rah=RaT.*Hall;
etar=eta(Ttop/Tb)./eta(Tint/Tb);
etar(find(isnan(etar)))=0;
% etat=exp(2.23);

if and(min(etar)<etat,max(etar)>etat);
    try
        bcv=interp1(etar-etat,ball,0);
    catch
        bcv=NaN;
    end
    Tiv=interp1(ball,Tint,bcv);
    if Tiv<=Tb; %no melting
        Ttv=interp1(ball,Ttop,bcv);
        DTv=Tb-Ttv;
        
        Hiv=H0*hfun(Tiv/Tb);
        Ftv=Fb+bcv.*Hiv;
        bbv=k0*(Tb-Tiv)/Fb;
        btv=k0*kfun(Ttv./Tb).*(Tiv-Ttv)./Ftv;
        biv=bcv-bbv-btv;
        if biv<0; %Not enough space in cell: no convection
            bcv=0;  bbv=0; btv=0; biv=0; 
            Tiv=Tb; Ttv=Tb; DTv=0; Hiv=H0; Ftv=Fb;etai=eta0;;etav=eta0;
        else
            etai=eta0*eta(Tiv/Tb);
            etav=eta0*eta(Ttv/Tb);
        end
    else
        bcv=0; bbv=0; btv=0; biv=0;
        Tiv=Tb; Ttv=Tb; DTv=0; Hiv=H0; Ftv=Fb;etai=eta0;;etav=eta0;
    end
else %never get the right viscosity
    bcv=0; bbv=0; btv=0; biv=0; 
    Tiv=Tb; Ttv=Tb; DTv=0; Hiv=H0; Ftv=Fb;etai=eta0;;etav=eta0;
end

