function [zout,tout,fout,kout,hout,d,r]=ConductionSolveSpread(dz,tb,fb,vn,kfun,hfun,ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [zout,tout,fout,kout,hout,d]=ConductionSolve(dz,tb,fb,vn,nz)
% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dz:   node spacing for FD grid
% tb:   bottom temperature
% fb:   bottom heat flux
% vn:   upwelling velocity
% kfun: Thermal conductivity @(t)
% hfun: heat production conductivity @(t/tb)
% ns:   number of nodes in desired output grid
% note: assumes that surface is at tb-1 (no-dimensional)
% output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zout: grid
% tout: temperature
% fout: heat flux
% kout: conductivity
% hout: heat production
% d:    shell thickness
% r:    integral of temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup analytical approximation for starting the solve
if vn==0
    ta=@(z,v,fb,tb)tb-fb*z-z.^2/2;
else
    ta=@(z,v,fb,tb)tb+z./v+(fb+1./v).*(1-exp(v.*z))./v;
end
% initialize first point
t=tb; k=kfun(1); h=hfun(1);
% ghost point
%tg=tb+fb*dz;
tg=ta(-dz,vn,fb,tb);
% move to second point
%t(2)=tb-fb*dz; 
t(2)=ta(dz,vn,fb,tb);
k(2)=kfun(t(2)/tb); h(2)=hfun(t(2)/tb);
% Handle third point
i=1;
t(i+2)=t(i)...
    +(2*dz*vn/k(i+1))*(t(i+1)-t(i))...
    +(k(i)/k(i+1))*(t(i+1)-tg)...
    -(h(i+1)+h(i))*dz^2/k(i+1);
k(i+2)=kfun(t(i+2)/tb); h(i+2)=hfun(t(i+2)/tb);

while ge(t(i+2),tb-1); 
    i=i+1;
    t(i+2)=t(i)...
        +(2*dz*vn/k(i+1))*(t(i+1)-t(i))...
        +(k(i)/k(i+1))*(t(i+1)-t(i-1))...
        -(h(i+1)+h(i))*dz^2/k(i+1);
    k(i+2)=kfun(t(i+2)/tb); h(i+2)=hfun(t(i+2)/tb);
end
% add one more point to allow flux calculation
i=i+1;
t(i+2)=t(i)...
    +(2*dz*vn/k(i+1))*(t(i+1)-t(i))...
    +(k(i)/k(i+1))*(t(i+1)-t(i-1))...
    -(h(i+1)+h(i))*dz^2/k(i+1);
%k(i+2)=kfun(t(i+2)/tb); h(i+2)=hfun(t(i+2)/tb);



nz=i+2;
f=[fb,-k(2:nz-1).*(t(3:nz)-t(1:nz-2))/(2*dz)];
z=0:nz-1;
d=interp1(t,z,tb-1)*dz;
r=integral(@(x)interp1(z*dz,t,x)-tb,0,d);

zout=linspace(0,d,ns);
tout=interp1(z,t,zout/dz);
fout=interp1(0:nz-2,f,zout/dz);

kout=tout*0+1;
hout=kout;

