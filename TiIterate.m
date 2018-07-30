function [DTn,Ran,Thn,Tcn,Tin]=TiIterate(bn,k0,K0,eta0,H0,Tb,Fb,kfun,eta,hfun,DT,Ra,Th)

%% start Ti=Tb, go on;
Tin=Tb; Ti=Tin; goon=1; it=0;
while goon; %abs(Ti-Tin)>0.01;
    Ti=Tin; it=it+1;
    kn=k0*kfun(Tin/Tb);
    Kn=K0*kfun(Tin/Tb);
    etan=eta0*eta(Tin/Tb);
    Hn=H0*hfun(Tin/Tb);
    
    DTn=DT(bn,kn,Kn,etan,Hn,Fb);
    Ran=Ra(bn,DTn,Kn,etan);
    Thn=Th(bn,DTn,Hn,kn,Ran);
    Tcn=Tb-DTn;
    Tin=Tb-(1-Thn)*DTn;
    if Tin<=0;
        Tin=0;DTn=Tb;Ran=NaN;Thn=NaN;Tcn=NaN;
        goon=0;
    end
    if abs(Ti-Tin)<0.01 
        goon=0;
    end
    if it>100;
        Tin=NaN;DTn=NaN;Ran=NaN;Thn=NaN;Tcn=NaN;
        goon=0;
    end
    % pn=eta(Tcn/Tb)/eta(Tin/Tb)
end

return