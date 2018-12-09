%% Preamble
% Each time starting new Matlab session execute in Matlab command window:
% irf
%%%%%%%%%%%%%%%%%%%%%%%To do: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
clear all
plt = 0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 for plotting
dt1sec  = 1.15740986075252e-05;
dt = 4*dt1sec;
global tfgm tedp tfpie tfpii tstate r B E ve vi Ne Ni Te Ti yyyymmddHHMM starttime endtime cavitystart cavityend; % variables from ParseData()
global c e eps0  kB mu0 me mp;
constants();
SetDefaultInterpretor();
instrument          = {'fgm','edp','fpi'};
[fgm,edp,fpii,fpie] = LoadData(instrument);
ParseData(fgm, edp, fpie, fpii);
%% Define time interval to perform Hybrid MVA
disp(yyyymmddHHMM)
formatout = 'dd-mmm-yyyy';
t1        = datenum([datestr(datenum(yyyymmddHHMM,'yyyymmdd'),formatout) ' ' starttime]);
t2        = datenum([datestr(datenum(yyyymmddHHMM,'yyyymmdd'),formatout) ' ' endtime]);
tcavity1  = datenum([datestr(datenum(yyyymmddHHMM,'yyyymmdd'),formatout) ' ' cavitystart]);
tcavity2  = datenum([datestr(datenum(yyyymmddHHMM,'yyyymmdd'),formatout) ' ' cavityend]);
[L,M,N]   = HybridMinVarAnalysis(fgm,B,t1,t2);
% rotate quantities in GSE xyz to LMN
[BLMN] = irf_newxyz(B(:,1:3), L, M, N); % B-field
[ELMN] = irf_newxyz(E(:,1:3), L, M, N); % E-field
[veLMN]= irf_newxyz(ve(:,1:3), L, M, N); % electron velocity
[viLMN]= irf_newxyz(vi(:,1:3), L, M, N); % ion velocity
%% Current density from FPI
% need to interpolate the ion velocity and ion density to match the
% number of data points as the electrons data. 
n = length(ve);
viLMNinterp(:,1)= interpolation(viLMN(:,1),n); 
viLMNinterp(:,2)= interpolation(viLMN(:,2),n);
viLMNinterp(:,3)= interpolation(viLMN(:,3),n);
Niinterp(:,1)   = interpolation(Ni,n);
nivi(:,1)       = Niinterp.*viLMNinterp(:,1);
nivi(:,2)       = Niinterp.*viLMNinterp(:,2);
nivi(:,3)       = Niinterp.*viLMNinterp(:,3);
neve(:,1)       = Ne.*veLMN(:,1);
neve(:,2)       = Ne.*veLMN(:,2);
neve(:,3)       = Ne.*veLMN(:,3);
JLMN(:,1)       = e*1e18*(nivi(:,1)-neve(:,1)); %1e18=1e9(nA)*1e3(m)*1e6(m^-3)
JLMN(:,2)       = e*1e18*(nivi(:,2)-neve(:,2));
JLMN(:,3)       = e*1e18*(nivi(:,3)-neve(:,3));
%% Current density from curlometer technique

%% Important indices in timeseries
% indices bracketing exhaust interval and cavity interval for fgm, fpi, edp instruments. 
% ix1 = start of exhaust, ix2 = end of exhaust, ix1cav = start of cavity,
% ix12cav = end of cavity
ix1i    = Index(fpii{1,1},t1);
ix2i    = Index(fpii{1,1},t2);
ix1icav = Index(fpii{1,1},tcavity1);
ix2icav = Index(fpii{1,1},tcavity2);
ix1e    = Index(fpie{1,1},t1);
ix2e    = Index(fpie{1,1},t2);
ix1ecav = Index(fpie{1,1},tcavity1);
ix2ecav = Index(fpie{1,1},tcavity2);
ix1B    = Index(fgm{1,1},t1);
ix2B    = Index(fgm{1,1},t2);
%% Check if jet is Alfvenic:
% select jet interval 
vijet   = viLMN(ix1i:ix2i,1:3);
% find size of jet (max v - average on either side of interval)
maxvi   = max(vijet(:,1)); % taking the max of +ve exhaust vL
minvi   = min(vijet(:,1));% taking the min of -ve exhaust vL
%maxvijet2= mean(vijet(:,1)); % taking the mean of the exhaust vL
viLin   = mean([mean(viLMN(ix1i-1:ix1i,1)),mean(viLMN(ix2i:ix2i+1,1))]);% use 2 data points left and right of interval
vijet   = max([abs(maxvi - viLin),abs(minvi - viLin)]); % magnitude of jet speed relative to inflow (take the bigger one)
%vijet1 = maxvijet2 - viLin; % speed of jet relative to inflow using mean
viin    = mean([mean(vi(ix1i-1:ix1i,4)),mean(vi(ix2i:ix2i+1,4))]);
% local Alfven speed- use average reconnecting field BL and |N| either side of the interval
BLin    = mean([abs(mean(BLMN(ix1B-1:ix1B,1))),abs(mean(BLMN(ix2B:ix2B+1,1)))]); %reconnection B-field (BR)
BMin    = mean([abs(mean(BLMN(ix1B-1:ix1B,2))),abs(mean(BLMN(ix2B:ix2B+1,2)))]); %Guide field (BG)
Niin    = mean([mean(Ni(ix1i-1:ix1i)),mean(Ni(ix2i:ix2i+1))]);
vAi     = BLin*10^(-9)/sqrt(mu0*Niin*10^6*mp)*10^(-3);
% ratio of the jet speed and Alfven speed 
jetratio = vijet/vAi; 
%% Walen Test

% correlation coefficient
%vprime = v - v_HT;
%CC = corrcoef(vprime, vAi); 
%% Heating - check the change in temperature in the exhaust relative to the inflow region 
ixTi = [ix1i,ix2i,ix1icav,ix2icav];
ixTe = [ix1e,ix2e,ix1ecav,ix2ecav];

% total TEMP of exhaust, cavity, and inflow regions
Tiexh = max(1/3*Ti(ixTi(1):ixTi(3),1) + 2/3*Ti(ixTi(1):ixTi(3),2));
Ticav = max(1/3*Ti(ixTi(3):ixTi(4),1) + 2/3*Ti(ixTi(3):ixTi(4),2));
Tiin1 = mean(1/3*Ti(ixTi(1)-1:ixTi(1),1) + 2/3*Ti(ixTi(1)-1:ixTi(1),2)); % use 2 data points left of interval
Tiin2 = mean(1/3*Ti(ixTi(2):ixTi(2)+1,1) + 2/3*Ti(ixTi(2):ixTi(2)+1,2)); % use 2 data points right of interval
Tiin  = mean([Tiin1,Tiin2]);
Teexh = mean(1/3*Te(ixTe(1):ixTe(3),1) + 2/3*Te(ixTe(1):ixTe(3),2));
Tecav = min(1/3*Te(ixTe(3):ixTe(4),1) + 2/3*Te(ixTe(3):ixTe(4),2));
Tein1 = mean(1/3*Te(ixTe(1)-1:ixTe(1),1) + 2/3*Te(ixTe(1)-1:ixTe(1),2)); % use 2 data points left of interval
Tein2 = mean(1/3*Te(ixTe(2):ixTe(2)+1,1) + 2/3*Te(ixTe(2):ixTe(2)+1,2)); % use 2 data points right of interval
Tein  = mean([Tein1,Tein2]);

% par and perp TEMP of exhaust, cavity, and inflow regions
Tiexhpar = max(Ti(ixTi(1):ixTi(3),1));
Ticavpar = max(Ti(ixTi(3):ixTi(4),1));
Tiexhperp = max(Ti(ixTi(1):ixTi(3),2));
Ticavperp = max(Ti(ixTi(3):ixTi(4),2));
TiparIn  = mean([mean(Ti(ixTi(1)-1:ixTi(1),1)),mean(Ti(ixTi(2):ixTi(2)+1,1))]);% use 2 data points left and right of interval
TiperpIn = mean([mean(Ti(ixTi(1)-1:ixTi(1),2)),mean(Ti(ixTi(2):ixTi(2)+1,2))]);% use 2 data points left and right of interval
Teexhpar = mean(Te(ixTe(1):ixTe(3),1));
Tecavpar = min(Te(ixTe(3):ixTe(4),1));
Teexhperp = mean(Te(ixTe(1):ixTe(3),2));
Tecavperp = min(Te(ixTe(3):ixTe(4),2));
TeparIn  = mean([mean(Te(ixTe(1)-1:ixTe(1),1)),mean(Te(ixTe(2):ixTe(2)+1,1))]);% use 2 data points left and right of interval
TeperpIn = mean([mean(Te(ixTe(1)-1:ixTe(1),2)),mean(Te(ixTe(2):ixTe(2)+1,2))]);% use 2 data points left and right of interval

% total heating for ion and electron
dTiexh = Tiexh - Tiin;  % speed of jet relative to inflow
dTicav = Ticav - Tiin;
dTiexhTheory = 0.13*mp*(vAi*10^3)^2/e; % (Drake et al., 2009; Phan et al., 2014)
dTeexh = Teexh - Tein;
dTecav = Tecav - Tein;
dTeexhTheory = 0.017*mp*(vAi*10^3)^2/e; % (Drake et  al., 2009; Phan et al., 2014)

% anisotropic heating
dTiexhpar = Tiexhpar-TiparIn;
dTicavpar = Ticavpar-TiparIn;
dTiexhperp = Tiexhperp-TiperpIn;
dTicavperp = Ticavperp-TiperpIn;
dTiparTheory = mp*(vAi*10^3)^2*BLin^2/(BLin^2+BMin^2)/e;
dTeexhpar = Teexhpar-TeparIn;
dTecavpar = Tecavpar-TeparIn;
dTeexhperp = Teexhperp-TeperpIn;
dTecavperp = Tecavperp-TeperpIn;

%% total plasma beta
Bin = mean([abs(mean(B(ix1B-1:ix1B,4))),abs(mean(B(ix2B:ix2B+1,4)))]);
beta = (Niin*10^6*(Tiin + Tein)*e)/((Bin*10^(-9))^2/(2*mu0));
betacrit = (viLin/vAi*pi*sqrt(2))^2; % assuming Zi=1 and mi = mp

%% Inflow conditions
Nein     = mean([mean(Ne(ix1e-1:ix1e)),mean(Ne(ix2e:ix2e+1))]);
BGoverBR = BMin/BLin;
B1       = mean(B(ix1B-1:ix1B,1:3));
B2       = mean(B(ix2B:ix2B+1,1:3));
Bshear   = acosd(dot(B1,B2)/(norm(B1)*norm(B2)));
wpi      = sqrt(Niin*10^6*e^2/mp/eps0);
wpe      = sqrt(Nein*10^6*e^2/me/eps0);
di       = c/wpi/10^3; % ion inertial length 
de       = c/wpe/10^3; % electron inertial length 
vgsein   = (mean(vi(ix1i-100:ix1i,1:3))+mean(vi(ix2i:ix2i+100,1:3)))/2; % average inflow vgse
exhdur   = (t2-t1)/dt1sec; % exhaust duration 
viNin    = abs(mean([mean(viLMN(ix1i-1:ix1i,3)),mean(viLMN(ix2i:ix2i+1,3))])); % average viN flow speed 
exhwidth = viNin*exhdur;
%% WRITE TO FILE
outputs = strcat(yyyymmddHHMM,',',starttime,',',endtime,',',cavitystart,',',cavityend,',',...
    num2str(Bin),',',num2str(Niin),',',num2str(Tiin),',',num2str(Tein),',',num2str(beta),',',num2str(betacrit),',',...
    num2str(BGoverBR),',',num2str(Bshear),',',num2str(L),',',num2str(M),',',num2str(N),',',num2str(vAi),',',num2str(vijet),',',num2str(vgsein),',',...
    num2str(dTiexh),',',num2str(dTiexhpar),',',num2str(dTiexhperp),',',num2str(dTiexhTheory),',',num2str(dTiparTheory),',',num2str(dTicav),',',num2str(dTicavpar),',',...
    num2str(dTicavperp),',',num2str(dTeexh),',',num2str(dTeexhpar),',',num2str(dTeexhperp),',',num2str(dTeexhTheory),',',num2str(dTecav),',',num2str(dTecavpar),',',...
    num2str(dTecavperp),',',num2str(di),',',num2str(de),',',num2str(exhwidth),',','\n');

fid = fopen('Results.csv','a+');
fprintf(fid,outputs);
fclose(fid);
% fprintf('vAi  v_in  vL_in  Ni_in  B_in  BL_in  BG/BR  B-shear  beta\n');
% fprintf('%4.1f %4.1f %5.1f %5.1f %5.1f %6.1f %6.1f %8.1f %5.1f\n', vAi, viin, viLin, Niin, Bin, BLin, BGoverBR, Bshear, beta)
%% Plotting 1
if plt==1
    make_it_tight = true;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.04], [0.1 0.01]);
    if ~make_it_tight,  clear subplot;  end

    % just plot +-4sec either side of exhaust interval 
    tfgm  = tfgm(Index(fgm{1,1},t1-dt):Index(fgm{1,1},t2+dt));
    tfpii = tfpii(Index(fpii{1,1},t1-dt):Index(fpii{1,1},t2+dt));
    tfpie = tfpie(Index(fpie{1,1},t1-dt):Index(fpie{1,1},t2+dt));
    tedp  = tedp(Index(edp{1,1},t1-dt):Index(edp{1,1},t2+dt));
    BLMN  = BLMN(Index(fgm{1,1},t1-dt):Index(fgm{1,1},t2+dt),1:3);
    viLMN = viLMN(Index(fpii{1,1},t1-dt):Index(fpii{1,1},t2+dt),1:3);
    veLMN = veLMN(Index(fpie{1,1},t1-dt):Index(fpie{1,1},t2+dt),1:3);
    Ni    = Ni(Index(fpii{1,1},t1-dt):Index(fpii{1,1},t2+dt));
    Ne    = Ne(Index(fpie{1,1},t1-dt):Index(fpie{1,1},t2+dt));
    Ti    = Ti(Index(fpii{1,1},t1-dt):Index(fpii{1,1},t2+dt),1:2);
    Te    = Te(Index(fpie{1,1},t1-dt):Index(fpie{1,1},t2+dt),1:2);
    ELMN  = ELMN(Index(edp{1,1},t1-dt):Index(edp{1,1},t2+dt),1:3);
    JLMN  = JLMN(Index(fpie{1,1},t1-dt):Index(fpie{1,1},t2+dt),1:3);

    % Init figure
    figure('Name','Replicate Jonathan plot');
    subplot(11,1,1);
    plot(tfgm,BLMN);
    ylabel({'B';'[nT]'},'Rotation',0);
    title('MMS3');
    legend('L','M','N');

    subplot(11,1,2);
    plot(tfpii,viLMN(:,1));
    hold on;
    plot(tfpie,veLMN(:,1));
    hold off;
    legend('viL','veL');
    ylabel({'$v_L$';'[km/s]'},'Rotation',0);

    subplot(11,1,3);
    plot(tfpii,viLMN(:,2));
    hold on;
    plot(tfpie,veLMN(:,2));
    hold off;
    legend('viM','veM');
    ylabel({'$v_M$';'[km/s]'},'Rotation',0);

    subplot(11,1,4);
    plot(tfpii,viLMN(:,3));
    hold on;
    plot(tfpie,veLMN(:,3));
    legend('viN','veN');
    ylabel({'$v_N$';'[km/s]'},'Rotation',0);

    subplot(11,1,5);
    plot(tfpii,Ni);
    hold on;
    plot(tfpie,Ne);
    legend('$N_i$','$N_e$');
    ylabel({'N';'[$cm^{-3}$]'},'Rotation',0);

    subplot(11,1,6);
    plot(tfpii,Ti(:,1));
    hold on;
    plot(tfpii,Ti(:,2));
    legend('Para','Perp');
    ylabel({'$T_i$';'[eV]'},'Rotation',0);

    subplot(11,1,7);
    plot(tfpie,Te(:,1));
    hold on;
    plot(tfpie,Te(:,2));
    legend('Para','Perp');
    ylabel({'$T_e$';'[eV]'},'Rotation',0);

    subplot(11,1,8);
    plot(tedp,ELMN);
    legend('L','M','N');
    ylabel({'E';'[mV/m]'},'Rotation',0);
    samexaxis();

    subplot(11,1,10);
    plot(tfpie,JLMN(:,1));
    legend('JL FPI');
    ylabel({'J';'[nA/$m^2$]'},'Rotation',0);
    samexaxis();

    %saveas(gcf,[yyyymmddHHMM '.png'])
    %saveas(gcf,[yyyymmddHHMM '.fig'])
end