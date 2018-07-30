%% parameters from Sotin and Labrosse
% Nu-Ra
c=0.3446;
xi=1/3;
zeta=4/3;
% internal temperature
e=1/2;
a=1.236;
beta=3/4;
gamma=1/4;
% c=(1e3)^(-xi)*e.^(-zeta)

% physical parameters
alpha=1.6e-4; %5.1e-5; % Expansivity [1/K]
g=1.315;        % Gravity [m^2/s]
eta0=10^(15);  % Viscosity prefactor [Pa/s]
rhoi=920;   % ice density   [kg/m^3]
rhow=1000;  % waterdensity  [kg/m^3]
k0=2.1;     % conductivity  [W/m/K]
K0=1.4e-6;  % Diffusivity   [m^2/s]
Tb=270;     % Basal temperature [K]
Ts=90;      % Surface temperature [K]
Fb=7e-3;    % Heat flux     [W/m^2]
Lh=333e3;   % Latent heat   [J/kg]
A0=26;      % scaled activation energy (from Sam)
edot=2e-10; % Strain rate   [1/s]
omega=2.048e-5; % Orbital frequency [1/s]
mu=4e9;     % Shear modulus [Pa]
Omega0=omega*eta0/mu;
H0=2*edot^2*eta0/(1+Omega0^2); % heat production at T=Tb, used for scaling [W/m^3]

kfun=@(tp)(1./tp);   % conductivity function (non-dim, consistent with k0 at melting
eta=@(tp)exp(A0.*(1./tp-1)); %t'=T/Tb viscosity dependence
hfun=@(tp)(1+Omega0^2).*eta(tp)./(1+(Omega0*eta(tp)).^2); %functional form

E0=100;     % Target band elevation [m]
Di=020e3;        % initial/maximum ice shell thickness to consider
nz=100; %number of points in export grid

Vall=10.^[-12:.05:-8]; %velocities to consider [m/s]
nV=numel(Vall);

%% Scales
T0=Tb-Ts;
L0=sqrt(k0*T0/H0);
V0=K0/L0;
F0=H0*L0;


rlv0=rhoi*Lh*V0;    % freezing

%% convection scaling relations
DT=@(b,k,K,eta,H,Fb)(1/e)*k^(-3/4)*(alpha*g*rhoi/K/eta)^(-1/4).*(((Fb+H*b)/c).^(3/4)-a*(H*b).^(3/4));
Ra=@(b,DT,K,eta)alpha*g*rhoi*DT.*b.^3/K/eta;
Th=@(b,DT,H,k,Ra)e+a*(H.*(b.^2)./k./DT).^beta.*(Ra.^(-gamma));
Hs=@(b,H,k,DT)H.*(b.^2)./k./DT;
Qt=@(b,H,Fc)Fb+H*b;%c.*(Ra(b).^xi).*(Th(b).^zeta)*k0.*DT(b)./b;

etat=exp(2.23);

Drft=NaN(1,nV);Vall;Dcrft=Drft;Tcrft=Drft;rrft=Drft;
Zrft=NaN(nz+3,nV);Trft=Zrft;Frft=Zrft;
% Zrft=NaN(nz,nV);Trft=Zrft;Frft=Zrft;
for iV=1:nV;
    vn=Vall(iV)/V0;
    Fn=Fb+vn*rlv0;
    
%     [zout,tout,fout,kout,hout,dn,rn]=ConductionSolveRift(1e-4,Tb/T0,Fn/F0,vn,kfun,hfun,100,2e4/L0);
    [zout,tout,fout,dn,rn,dcon,tcon]=solvestagnantspread(Vall(iV),Di,etat,k0,K0,eta0,H0,Tb,Ts,Fn,kfun,eta,hfun,Hs,DT,Ra,Th,T0,L0,F0); %%
%     [zout,tout,fout,dn,rn,dcon,tcon]=solvestagnantrift(vn,Di,etat,k0,K0,eta0,H0,Tb,Ts,Fn,kfun,eta,hfun,Hs,DT,Ra,Th,T0,L0,F0); %%
    % store dimensional
    Drft(iV)=dn*L0;
%     Dcrft(iV)=dn*L0;    Tcrft(iV)=Tb;
    Dcrft(iV)=dcon*L0;
    Tcrft(iV)=tcon*T0;
    Zrft(:,iV)=(zout)*L0; %also convert to depth
    Trft(:,iV)=tout*T0;
    Frft(:,iV)=fout*F0;
    rrft(iV)=rn;
    
%     [zout,tout,fout,kout,hout,dn,rn]=ConductionSolveRift(1e-4,Tb/T0,Fn/F0,vn,kfun,hfun,100,2e4/L0);
    [zout,tout,fout,kout,hout,dn,rn]=ConductionSolveSpread(1e-4,Tb/T0,Fn/F0,vn,kfun,hfun,100);
    Dcd(iV)=dn*L0;
    Zcd(:,iV)=(zout)*L0; %Also convert to depth
    Tcd(:,iV)=tout*T0;
    Fcd(:,iV)=fout*F0;
    rcd(iV)=rn;
end
%%
[znv,tnv,fnv,dnv,rnv,dcnv,tcnv]=solvestagnantstatic(Di,etat,k0,K0,eta0,H0,Tb,Ts,Fb,kfun,eta,hfun,Hs,DT,Ra,Th,T0,L0,F0); %% 
Dnv=dnv*L0;
Znv=znv*L0;
Tnv=tnv*T0;
Fnv=fnv*F0;
Dcnv=dcnv*L0;
Tcnv=tcnv*T0;

% %%
% figure(1); clf; hold on;
% 
% plot(tout*T0,zout*L0/1000,'linewidth',2);
% plot(tcon*T0,dcon*L0/1000,'o')
% % plot(Tnv,Znv/1000,'linewidth',2);
% 
% set(gca,'ydir','reverse','box','on','fontsize',12); 
% xlabel('Temperature (K)','fontsize',18);
% ylabel('Depth (km)','fontsize',18);
% 

%% elevation elements
% EthickSpr=(Dnv-Dspr).*(rhoi-rhow)./rhow;
% EthermSpr=(rnv-rspr)*alpha*T0*L0*rhoi/rhow;
% % required old ice density
% rho1Spr=(rhow*(E0+Dnv-Dspr)+rhoi*(Dspr+alpha*T0*L0*rspr))./(Dnv+alpha*T0*L0*rnv);
% rho2Spr=(rhoi*(E0+Dnv-Dspr)+rhoi*(Dspr+alpha*T0*L0*rspr))./(Dnv+alpha*T0*L0*rnv);
% rho3Spr=rhoi*(Dnv+E0)./Dnv;

EthickRft=(Dnv-Drft).*(rhoi-rhow)./rhow;
EthermRft=(rnv-rrft)*alpha*T0*L0*rhoi/rhow;
% EthickCd=(Dnv-Dcd).*(rhoi-rhow)./rhow;
% EthermCd=(rnv-rcd)*alpha*T0*L0*rhoi/rhow;
% required old ice density
rho1Rft=(rhow*(E0+Dnv-Drft)+rhoi*(Drft+alpha*T0*L0*rrft))./(Dnv+alpha*T0*L0*rnv);
rho2Rft=(rhoi*(E0+Dnv-Drft)+rhoi*(Drft+alpha*T0*L0*rrft))./(Dnv+alpha*T0*L0*rnv);
rho3Rft=rhoi*(Dnv+E0)./Dnv;
%%
figure(1); clf; 
subplot 211; hold on
% semilogx(Vall,[Drft;Dcrft;Drft-Dcrft]/1000,'LineWidth',2);
semilogx(Vall,[Drft;Drft-Dcrft]/1000,'LineWidth',2);
semilogx(min(Vall)*[1,1],[Dnv;Dnv-Dcnv]/1000,'ko')

% semilogx(Vall,[Drft;Drft-Dcrft;Dcd]/1000,'LineWidth',2);
set(gca,'fontsize',12,'box','on','xscale','log')
xlabel('Upwelling velocity (m/s)','fontsize',18);
ylabel('Thickness (km)','fontsize',18);
% legend('Ice thickness','Conductive ice thickness','Convective cell thickness','location','best');
legend('Ice thickness','Convective cell thickness','location','best');
% legend('Ice thickness','Convective cell thickness','Conduction only','location','best');

subplot 212
Fba=Fb+Vall*rlv0/V0;
Fca=Fba;
for iV=1:nV;
    if Dcrft(iV)<Drft(iV);
        Fca(iV)=interp1(Zrft(:,iV),Frft(:,iV),Dcrft(iV));
    end       
end
loglog(Vall,[Fba;Fca;Frft(end,:)]*1000,'LineWidth',2)
set(gca,'fontsize',12,'box','on')
set(gca,'YScale','linear','ylim',[0,50],'xlim',[min(Vall),max(Vall)])
xlabel('Upwelling velocity (m/s)','fontsize',18);
ylabel('Heat flux (mW/m^2)','fontsize',18);
legend('Base of the ice','Base of conducting ice','Surface','location','best');
%%
figure(2); clf; 
subplot(221); hold on;
Hr=plot(Trft,Zrft/1e3,'LineWidth',2); 
Hc=plot(Tcrft,Dcrft/1e3,'sk');

% Hs=plot(Tspr,Zspr/1e3,'--','LineWidth',2); 
% for i=1:numel(Hs);set(Hs(i),'color',get(Hr(i),'color'));end;
set(gca,'ydir','reverse','fontsize',12,'box','on')
xlabel('Temperature (K)','fontsize',18);
ylabel('Depth (km)','fontsize',18);
ylim([0,40]);

subplot(222); hold on;
Hr=plot(Frft*1000,Zrft/1e3,'LineWidth',2); 
% Hs=plot(Fspr*1000,Zspr/1e3,'--','LineWidth',2); 
% for i=1:numel(Hs);set(Hs(i),'color',get(Hr(i),'color'));end;
set(gca,'ydir','reverse','fontsize',12,'xscale','log','box','on')
xlabel('Heat flux (mW/m^2)','fontsize',18);
ylabel('Depth (km)','fontsize',18);
subplot 221
plot(Tnv,Znv/1e3,'k','linewidth',2)
subplot 222
plot(Fnv*1000,Znv/1e3,'k','linewidth',2)
if nV<=8;
    Hl=legend(num2str(Vall'),'location','SE');
    title(Hl,'Upwelling velocity');
end
ylim([0,40]);

subplot 413; hold on
Hr=plot(Vall,[EthickRft+EthermRft],'k','linewidth',2);% (rhoi/rhow-1)*(Dnv-Dspr));
% Hrcd=plot(Vall,[EthickCd+EthermCd],'b--','linewidth',2);% (rhoi/rhow-1)*(Dnv-Dspr));
% Hs=plot(Vall,[EthickSpr+EthermSpr],'k--','linewidth',2);% (rhoi/rhow-1)*(Dnv-Dspr));
% for i=1:numel(Hs);set(Hs(i),'color',get(Hr(i),'color'));end;
% legend('Rifting','Spreading','location','SW')
set(gca,'XScale','log','box','on','FontSize',12)
xlabel('Velocity (m/s)','fontsize',18);
ylabel('Elevation (m)','fontsize',18);
legend(sprintf('Maximum elevation: %g m',max ([EthickRft+EthermRft])))
% legend(sprintf('Maximum elevation with convection: %g m',max ([EthickRft+EthermRft])),...
%     sprintf('Maximum elevation without convection: %g m',max ([EthickCd+EthermCd])))

subplot 414; hold on;
% hb=100;
%plot(Vall,rhow*(1+hb/Dnv)+(rhoi-rhow)*Dspr/Dnv);
Hr=plot(Vall,[rho1Rft;rho2Rft;rho3Rft*ones(size(Vall))],'LineWidth',2); 
plot([min(Vall),max(Vall)],1000*[1,1],'k','linewidth',0.5)
% Hs=plot(Vall,[rho1Spr;rho2Spr;rho3Spr*ones(size(Vall))],'--','LineWidth',2); 
% for i=1:numel(Hs);set(Hs(i),'color',get(Hr(i),'color'));end;
legend('\rho_1','\rho_2','\rho3','location','NW')
set(gca,'XScale','log','box','on','FontSize',12)
ylim([920,1020])
xlabel('Velocity (m/s)','fontsize',18);
ylabel('Ice density (kg/m^3)','fontsize',18);

return
%%
% Tint=Tb-(1-Th(ball)).*(DT(ball));
% Ttop=Tb-DT(ball);
% Tbot=Tb+ball*0;
figure(2);clf;
subplot 221; hold on;
%plot(ball./1000, DT(ball),'linewidth',2);
% set(gca,'box','on','fontsize',12);
% xlabel('Convective thickness [km]','fontsize',18);
% ylabel('Temperature drop [K])','fontsize',18);
plot(ball./1000, Tint,'linewidth',2);
plot(ball./1000, Ttop,'linewidth',2);
plot(ball./1000, Tbot,'--','linewidth',2)
hl=legend('Interior','Top','Bottom','location','SW'); set(hl,'fontsize',18)
set(gca,'box','on','fontsize',12);
xlabel('Convective thickness [km]','fontsize',18);
%ylabel('Interior temperature [K])','fontsize',18);
ylabel('Temperature [K])','fontsize',18);

subplot 222; hold on;
plot(ball./1000, RaT,'linewidth',2);
% plot(ball./1000, Rah,'linewidth',2);
% plot(ball./1000, Hall,'linewidth',2);
YL=ylim; ylim([1,max(YL)]);

set(gca,'box','on','fontsize',12,'yscale','log');
xlabel('Convective thickness [km]','fontsize',18);
ylabel('Rayleigh Number','fontsize',18);
% ylabel('Number','fontsize',18);
% hl=legend('Ra_T','Ra_H','H_s','location','SE'); set(hl,'fontsize',18)

subplot 223; hold on;
%plot(ball./1000, Th(ball),'linewidth',2);
% plot(ball./1000, Tb-(1-Th(ball)).*(DT(ball)),'linewidth',2);
plot(ball./1000,eta0*eta(Tint./Tbot),'linewidth',2);
plot(ball./1000,eta0*eta(Ttop./Tbot),'linewidth',2);
plot(ball./1000,eta0*eta(Tbot./Tbot),'--','linewidth',2)

set(gca,'box','on','fontsize',12,'yscale','log');
xlabel('Convective thickness [km]','fontsize',18);
ylabel('Viscosity [Pa.s]','fontsize',18);
% ylim(eta0*10.^[-1,2])
% ylabel('Interior temperature [K])','fontsize',18);
hl=legend('Interior','Top','Bottom','location','NW'); set(hl,'fontsize',18)
subplot 224; hold on;
% %plot(ball./1000, Th(ball),'linewidth',2);
plot(ball./1000, (Fb+ball.*H0.*hfun(Tint/Tb)).*1000,'linewidth',2,'color',[0,0,0]);
plot([min(ball/1000),max(ball/1000)],Fb*[1,1]*1000','--','linewidth',2,'color',[0,0,0])
% plot(ball./1000, Qtall*1000,'--c','linewidth',2);

set(gca,'box','on','fontsize',12);
xlabel('Convective thickness [km]','fontsize',18);
ylabel('Heat flux [mW/m^2]','fontsize',18);
hl=legend('Top','Bottom','location','NW'); set(hl,'fontsize',18)
%%

etat=exp(2.23);
etar(find(isnan(etar)))=0;
if and(min(etar)<etat,max(etar)>etat);
    bcv=interp1(etar-etat,ball,0)
else
    bcv=0;
end
%
subplot 221
plot([0,1,1,0]*bcv/1000,interp1(ball,Ttop,bcv)*[0,0,1,1]+interp1(ball,Tint,bcv)*[1,1,0,0],'k')
hl=legend('Interior','Top','Bottom','location','SW'); set(hl,'fontsize',18)
subplot 222
plot([0,1,1]*bcv/1000,interp1(ball,RaT,bcv)*[1,1,0]+[0,0,1],'k')

% hl=legend('Ra_T','Ra_H','H_s','location','SE'); set(hl,'fontsize',18)
subplot 223
plot([0,1,1,0]*bcv/1000,eta0*eta(interp1(ball,Ttop,bcv)/Tb)*[0,0,1,1]+eta0*eta(interp1(ball,Tint,bcv)/Tb)*[1,1,0,0],'k')

hl=legend('Interior','Top','Bottom','location','NW'); set(hl,'fontsize',18)
subplot 224
plot([0,1,1]*bcv/1000,(interp1(ball,Qtall,bcv)*[1,1,0]+Fb*[0,0,1])*1000,'k')
hl=legend('Top','Bottom','location','NW'); set(hl,'fontsize',18)

%%

bb=k0*(Tbot-Tint)/Fb;
bt=k0*kfun(Ttop./Tbot).*(Tint-Ttop)./Qtall;
bi=ball-bb-bt;



return
%%
figure(3); clf; hold on;
G=2.23/A0
DTall=DT(ball);
Thall=Th(ball);
xall=(1-Thall).*DTall./Tb;
gall=Thall./G./(1-Thall);
DTp=(2+gall).*(1+sqrt(gamma./(2+gamma)));
DTm=(2+gall).*(1-sqrt(gamma./(2+gamma)));
plot(ball./1000, DTall,'linewidth',2);
plot(ball./1000, DTm,'linewidth',2);
plot(ball./1000, DTp,'linewidth',2);
ylim([0,Tb-Ts])



%% Target viscosity contrast:
etat=1000;
DTt=Tb*(1-1./(1+log(etat)./A0));
bvisc=fzero(@(b)DT(b)-DTt,[0,Di]);
bthi=fzero(@(b)Th(b)-1,[1,Di]);
b=min(bthi,bvisc)
Racv=Ra(b);
Qtcv=Qt(b);

%% Scales
T0=Tb-Ts;
L0=sqrt(k0*T0/H0);
V0=K0/L0;
F0=H0*L0;

% initialize
% Dspr=NaN(1,nV); Drft=Dspr;     % ice thickness
% Tspr=NaN(nz,nV); Trft=Tspr;    % temperature
% Fspr=NaN(nz,nV); Frft=Fspr;    % heat flux
% Zspr=NaN(nz,nV); Zrft=Zspr;    % depth
% rspr=NaN(1,nV); rrft=rspr;     % integrated temperature (ND)
tb=Tb/T0;           % bottom temperature (ND)
rlv0=rhoi*Lh*V0;    % freezing

% conductive part
dz=1e-4;
% add no upwelling solution
fb=Fb/F0;
[znv,tnv,fnv,kout,hout,dnv,rnv]=...
        ConductionSolve(dz,tb,fb,0,kfun,hfun,nz);
Dnv=dnv*L0;
Znv=(dnv-znv)*L0;
Tnv=tnv*T0;
Fnv=fnv*F0;
% with convection
[zcv,tcv,fcv,kout,hout,dcv,rcv]=...
        ConductionSolve(dz,tb-DTt/T0,Qtcv/F0,0,kfun,hfun,nz);
Dcv=dcv*L0;
Zcv=(dcv-zcv)*L0;
Tcv=tcv*T0;
Fcv=fcv*F0;
% %%
% b=1e4; % 10 km thick ice
% DT=100; %100 T drop
% Ra=(alpha*g*rhoi*DT*b^3/K0/eta0)
% h=H0*b^2/k0/DT
% Th=e+a*h^beta*Ra^(-gamma)
% q=c*Ra^xi*Th^zeta
% Qt=q*k0*DT/b
% Qb=Qt-H0*b

