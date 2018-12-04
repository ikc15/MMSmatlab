function ParseData(fgm,edp,fpie,fpii)

global tfgm tedp tfpie tfpii tstate r B E ve vi Ne Ni Ti Te;

% convert epoch time (matlab autoconvert to datenum) to datetime:
tfgm = datetime(fgm{1,1},'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');
tedp = datetime(edp{1,1},'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');
tfpie = datetime(fpie{1,1},'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');
tfpii = datetime(fpii{1,1},'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS');

% convert epoch time for the state of the spacecraft to datetime:
tstate = datetime(fgm{1,7},'ConvertFrom','datenum','TicksPerSecond',1000,'Format','dd-MMM-yyyy HH:mm:ss.SSS');

% Spacecraft Position 
r(:,1:4) = fgm{1,8}; %x, y, z, r (GSE- in m)

% Magnetic field (B) 
B(:,1:4) = fgm{1,2}; %Bx,By,Bz,B (GSE- in nT)

% Electric field (E)
E(:,1:3) = edp{1,2}; %Ex, Ey, Ez (GSE- in mV/m)

% Velocity (v) 
ve(:,1:4) = irf_abs(fpie{1,29}); %vex,vey,vez,ve (GSE- in km/s)
vi(:,1:4) = irf_abs(fpii{1,25}); %vix,viy,viz,vi (GSE- in km/s)

% Number density (N) 
Ni(:,1) = fpii{1,19}; %cm^-3
Ne(:,1) = fpie{1,23}; %cm^-3

% Temperature (T) 
Ti(:,1)=fpii{1,40}; %1 = para (in eV)
Ti(:,2)=fpii{1,41}; %2 = perp (in eV)
Te(:,1)=fpie{1,45}; %1 = para (in eV)
Te(:,2)=fpie{1,46}; %2 = perp (in eV)
end