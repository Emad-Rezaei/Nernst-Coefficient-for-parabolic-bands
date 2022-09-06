%This code is to compute Nernst coefficient for multi isotropic bands
%Developed by Emad
clear all
clc
disp('This code is to compute Nernst coefficient within a two-parabolic-band model')
k=8.617*10^-5; %Boltzman Constant in ev
q=-1.6*10^-19; %electron charge in c not ev
%For TM and TE vol is not taken into account
vol=304.5466; %unit cell volume of GaAs in au^3
conv=0.529177; %conversion from au^3 to Ang^3
vol=vol*conv^3*10^-30; %-30 vol in m3,-24 in cm3
theta=75; %angel between crystal c-axis and perpendicualr z-axis 
mupt=1000; %number of points between the initial and final chemical potential 
Ept=1000; %number of points between the initial and final energy
Hvec=[0.00,0.00,1.00]; %magnetic field vector
H=Hvec(3);%magnetic field along z-axis
m0=9.109*10^-31; %electron mass
T=300.00; %Temperature
tau0=10.00*10^-15; %relaxation time Si
re=0.00;
rpl=-0.00;
rph=-0.0;
h=6.626*10^-34 ; %planck's constant in j.s
hbar=h/(2*pi); %hbar
mc=0.02050*m0; % 0.130
mv=0.025*m0; %heavy hole 0.19
mhs=0.031*m0; %0.0326*m0
mvl=0.07*m0; %1.26
Ev=10.07; %bottom of conduction band GaAs WAN 6.8116 6.7563552
Evs=10.0018;
Evl=10.07;
Eg=1.563737; %1.9190461
Ec=Ev+Eg; %top of valance band
muf=0.00; %intrinsic chemical potential
%Initiation of variables
ne=zeros(mupt,1); %electron concentration
np=zeros(mupt,1); %hole concentration
nt=zeros(mupt,1); %total carrier concentration
mu=zeros(mupt,1); %chemical potential 
xf=zeros(mupt,1); %reduced chemical potential 
sigxx=zeros(mupt,1); %conduction band electrical conductivity  in xx-direction
kcoxx=zeros(mupt,1); %total electronic thermal conductivity 
kcoxyz=zeros(mupt,1); %total electronic thermal conductivity 2nd order xyz component
kcoyxz=zeros(mupt,1);
Bxx=zeros(mupt,1); %1st order B response function holes
Bxyz=zeros(mupt,1); %2nd order B response function conduction band xyz
Byxz=zeros(mupt,1);
sigyxz=zeros(mupt,1); %electrical conductivity holes YXZ
sigxyz=zeros(mupt,1); %total 2nd order electrical conductivity 
detsig=zeros(mupt,1); %determinant of total sigma
ITN=zeros(mupt,1); %Isothermal Nernst
ETN=zeros(mupt,1); %Ettinghausen
msig=zeros(3);
mB=zeros(3);
mkco=zeros(3,3);
mkap=zeros(mupt,3,3);
mrho=zeros(mupt,3,3);
mS=zeros(mupt,3,3); %Seebeck Energy xx yy zz
E=zeros(Ept,1); %Energy in eV
dosc=zeros(Ept,1); %DOS of conduction band
dosv=zeros(Ept,1); %DOS of valence band
tdfc=zeros(Ept,1); %1st order Transport distribution function conduction band
tdfvh=zeros(Ept,1); %1st order Transport distribution function valence band heavy hole
tdfvs=zeros(Ept,1); %1st order Transport distribution function valence band heavy hole up
tdfvl=zeros(Ept,1); %1st order Transport distribution function
TDF1=zeros(Ept,1); %total 1st order Transport distribution function
tdf2c=zeros(Ept,3); %2nd Transport distribution function
tdf2vh=zeros(Ept,3); %2nd Transport distribution function
tdf2vs=zeros(Ept,3); %2nd Transport distribution function valence up
tdf2vl=zeros(Ept,3); %2nd Transport distribution function
TDF2=zeros(Ept,3); %2nd Transport distribution function
f=zeros(Ept,mupt); %F-D distribution funtion
mindfdE=zeros(Ept,mupt); %minus derivative of F-D distribution funtion
Ac=inv(3*pi^2*hbar^3)*2*sqrt(2*mc)*tau0*q^2; %coefficient for 1st order TDF of conductoin band
Avh=inv(3*pi^2*hbar^3)*2*sqrt(2*mv)*tau0*q^2; %coefficient for 1st order TDF of valence band heavy
Avhs=inv(3*pi^2*hbar^3)*2*sqrt(2*mhs)*tau0*q^2; %coefficient for 1st order TDF of valence band heavy up
Avl=inv(3*pi^2*hbar^3)*2*sqrt(2*mvl)*tau0*q^2; %coefficient for 1st order TDF of valence band light
Ac2=0;
Avh2=0;
Avhs2=0;
Avl2=0;
if mc~=0
    Ac2=inv(3*hbar^3*pi^2*sqrt(mc))*2*(q)^3*H*tau0^2*sqrt(2); 
end
%- sing is required due to the nagative dv/dk at valence bands
if mv~=0
    Avh2=-inv(3*hbar^3*pi^2*sqrt(mv))*2*(q)^3*H*tau0^2*sqrt(2); %coefficient for 2nd order TDF of valence band 
end
if mhs~=0
    Avhs2=-inv(3*hbar^3*pi^2*sqrt(mhs))*2*(q)^3*H*tau0^2*sqrt(2); %coefficient for 2nd order TDF of valence band up
end
if mvl~=0
    Avl2=-inv(3*hbar^3*pi^2*sqrt(mvl))*2*(q)^3*H*tau0^2*sqrt(2); %coefficient for 2nd order TDF of valence band
end
%Energy grid generation
for i=1:size(E,1)
    E(i,1)=1/size(E,1)*i*(Ec+Ev);
    tdf2vh(i,1)=E(i,1); %energy
    tdf2vs(i,1)=E(i,1); %energy
    tdf2vl(i,1)=E(i,1); %energy
    tdf2c(i,1)=E(i,1); %energy
    TDF2(i,1)=E(i,1); %energy
    if E(i,1)<Ev
       tdfvh(i,1)=Avh*inv(k*T)^rph*(Ev-E(i,1))^rph*(Ev-E(i,1))^(1.5)*abs(q)^(1.5);%why abs(q)^(1.5)?E^1.5 in TDF which is in ev
       tdfvs(i,1)=Avhs*inv(k*T)^rph*(Ev-E(i,1))^rph*(Ev-E(i,1))^(1.5)*abs(q)^(1.5);%why abs(q)^(1.5)?E^1.5 in TDF which is in ev
       dosv(i,1)=mv*sqrt(2*mv*(Ev-E(i,1)))/(pi^2*hbar^3)*abs(q)^1.5; %q^1.5 conversion from eV to J 
       tdf2vh(i,2)=Avh2*inv(k*T)^(2*rph)*(Ev-E(i,1))^(2*rph)*(Ev-E(i,1))^1.5*abs(q)^(1.5); %TDFxyz because in smith -\Omega tau(E) not a constant (Ev-E(i,1))^r
       tdf2vh(i,3)=-Avh2*inv(k*T)^(2*rph)*(Ev-E(i,1))^(2*rph)*(Ev-E(i,1))^1.5*abs(q)^(1.5); %TDFyxz because in smith -\Omega tau(E) not a constant (Ev-E(i,1))^r
    end
    if E(i,1)<Evs
        tdf2vs(i,2)=Avhs2*inv(k*T)^(2*rph)*(Evs-E(i,1))^(2*rph)*(Evs-E(i,1))^1.5*abs(q)^(1.5);
        tdf2vs(i,3)=-Avhs2*inv(k*T)^(2*rph)*(Evs-E(i,1))^(2*rph)*(Evs-E(i,1))^1.5*abs(q)^(1.5);
    end
    if E(i,1)<Evl
       tdf2vl(i,2)=Avl2*inv(k*T)^(2*rpl)*(Evl-E(i,1))^(2*rpl)*(Evl-E(i,1))^1.5*abs(q)^(1.5); %TDFxyz because in smith -\Omega tau(E) not a constant (Ev-E(i,1))^r
       tdf2vl(i,3)=-Avl2*inv(k*T)^(2*rpl)*(Evl-E(i,1))^(2*rpl)*(Evl-E(i,1))^1.5*abs(q)^(1.5); %TDFyxz because in smith -\Omega tau(E) not a constant (Ev-E(i,1))^r
       tdfvl(i,1)=Avl*inv(k*T)^rpl*(Evl-E(i,1))^rpl*(Ev-E(i,1))^(1.5)*abs(q)^(1.5);%why abs(q)^(1.5)?E^1.5 in TDF which is in ev
    end
    if E(i,1)>Ec
       tdfc(i,1)=Ac*inv(k*T)^re*(E(i,1)-Ec)^re*(E(i,1)-Ec)^(1.5)*abs(q)^(1.5); %tau(E) not a constant (E(i,1)-Ec)^r
       dosc(i,1)=mc*sqrt(2*mc*(E(i,1)-Ec))/(pi^2*hbar^3)*abs(q)^1.5; %q^1.5 conversion from eV to J
       tdf2c(i,2)=Ac2*inv(k*T)^(2*re)*(E(i,1)-Ec)^(2*re)*(E(i,1)-Ec)^1.5*abs(q)^(1.5); %TDFxyz because in smith -\Omega tau(E) not a constant (E(i,1)-Ec)^r
       tdf2c(i,3)=-Ac2*inv(k*T)^(2*re)*(E(i,1)-Ec)^(2*re)*(E(i,1)-Ec)^1.5*abs(q)^(1.5); %TDFyxz because in smith -\Omega tau(E) not a constant (E(i,1)-Ec)^r
    end
    TDF2(i,2)=tdf2c(i,2)+tdf2vh(i,2)+tdf2vs(i,2)+tdf2vl(i,2);
    TDF2(i,3)=tdf2c(i,3)+tdf2vh(i,3)+tdf2vs(i,3)+tdf2vl(i,3);
end
TDF2(:,2)=tdf2c(:,2)+tdf2vh(:,2)+tdf2vs(:,2)+tdf2vl(:,2);
TDF2(:,3)=tdf2c(:,3)+tdf2vh(:,3)+tdf2vs(:,3)+tdf2vl(:,3);
TDF1=tdfvh+tdfvs+tdfvl+tdfc;
%loop over chemical potential
for j=1:size(mu,1)
    mu(j,1)=0.5*(Ec+Ev)-1+2*j/size(mu,1); 
    xf(j,1)=inv(k*T)*(mu(j,1)-0.5*Eg-Ev);
    for i=1:size(E,1)
        if E(i,1)<=mu(j,1)+k*T*(log(inv(1e-9)-1))
            f(i,j)=inv(1+exp((E(i,1)-mu(j,1))/(k*T)));
        end
        mindfdE(i,j)= exp((E(i,1)-mu(j,1))/(k*T))/(k*T*(1+exp((E(i,1)-mu(j,1))/(k*T)))^2); 
    end
end
for j=1:size(mu,1)
    for i=1:size(E,1)-1
        if E(i,1)>=Ec
            ne(j,1)=ne(j,1)+(dosc(i,1)*f(i,j)+dosc(i+1,1)*f(i+1,j))*0.5*(E(i+1,1)-E(i,1));
        end
        if E(i,1)<=Ev
            np(j,1)=np(j,1)+(dosv(i,1)*(1-f(i,j))+dosv(i+1,1)*(1-f(i+1,j)))*0.5*(E(i+1,1)-E(i,1));
        end
        if abs(ne(j,1)-np(j,1))<=1e-3
            muf=mu(j,1);
        end
        %trpsignxx(j,1)=trapz(E,tdfc.*mdfdE(:,j));
        sigxx(j,1)=sigxx(j,1)+0.5*(E(i+1,1)-E(i,1))*(TDF1(i+1,1)*mindfdE(i+1,j)+TDF1(i,1)*mindfdE(i,j));
        Bxx(j,1)=Bxx(j,1)+0.5*inv(-T)*(E(i+1,1)-E(i,1))*(TDF1(i+1,1)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))+TDF1(i,1)*mindfdE(i,j)*(E(i,1)-mu(j,1)));
        kcoxx(j,1)=kcoxx(j,1)+0.5*inv(T)*(E(i+1,1)-E(i,1))*(TDF1(i+1,1)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))^2+TDF1(i,1)*mindfdE(i,j)*(E(i,1)-mu(j,1))^2);
        sigxyz(j,1)=sigxyz(j,1)+0.5*(E(i+1,1)-E(i,1))*(TDF2(i+1,2)*mindfdE(i+1,j)+TDF2(i,2)*mindfdE(i,j));
        Bxyz(j,1)=Bxyz(j,1)+0.5*inv(-T)*(E(i+1,1)-E(i,1))*(TDF2(i+1,2)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))+TDF2(i,2)*mindfdE(i,j)*(E(i,1)-mu(j,1)));
        kcoxyz(j,1)=kcoxyz(j,1)+0.5*inv(T)*(E(i+1,1)-E(i,1))*(TDF2(i+1,2)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))^2+TDF2(i,2)*mindfdE(i,j)*(E(i,1)-mu(j,1))^2);
        sigyxz(j,1)=sigyxz(j,1)+0.5*(E(i+1,1)-E(i,1))*(TDF2(i+1,3)*mindfdE(i+1,j)+TDF2(i,3)*mindfdE(i,j));
        Byxz(j,1)=Byxz(j,1)+0.5*inv(-T)*(E(i+1,1)-E(i,1))*(TDF2(i+1,3)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))+TDF2(i,3)*mindfdE(i,j)*(E(i,1)-mu(j,1)));
        kcoyxz(j,1)=kcoyxz(j,1)+0.5*inv(T)*(E(i+1,1)-E(i,1))*(TDF2(i+1,3)*mindfdE(i+1,j)*(E(i+1,1)-mu(j,1))^2+TDF2(i,3)*mindfdE(i,j)*(E(i,1)-mu(j,1))^2);
    end
    msig(1,1)=sigxx(j,1);
    msig(2,2)=sigxx(j,1);
    msig(3,3)=sigxx(j,1);
    msig(1,2)=sigxyz(j,1);
    msig(2,1)=sigyxz(j,1);
    mB(1,1)=Bxx(j,1);
    mB(2,2)=Bxx(j,1);
    mB(3,3)=Bxx(j,1);
    mB(1,2)=Bxyz(j,1);
    mB(2,1)=Byxz(j,1);
    mkco(1,1)=kcoxx(j,1);
    mkco(2,2)=kcoxx(j,1);
    mkco(3,3)=kcoxx(j,1);
    mkco(1,2)=kcoxyz(j,1);
    mkco(2,1)=kcoyxz(j,1);
    mS(j,:,:)=inv(msig)*mB;
    mrho(j,:,:)=inv(msig);
    mkap(j,:,:)=mkco-msig*((inv(msig)*mB)*(inv(msig)*mB))*T;
    ITN(j,1)=mS(j,2,1);
    %ITN(j,2)=(sig(j,2)*Byxz(j,2)-sigyxz(j,2)*B(j,2))/(sig(j,2)*sig(j,3)-sigxyz(j,2)*sigyxz(j,2));
    ITH(j,1)=mrho(j,2,1);
    ETN(j,1)=mS(j,2,1)*T/mkap(j,2,2);
end
%semilogy(mu(:,1),sigxx(:,1))
plot(mu(:,1)-Ev,ITN(:,1)*1e6*1.5)
%plot(mu(:,1)-0.5*(Ec+Ev),mS(:,1,1)*1e6)
%plot(mu(:,1),sigxx(:,1))
%xlim([-1.0 1.0]);
 %ne=ne/vol;
 %np=np/vol;
 nt=ne-np;