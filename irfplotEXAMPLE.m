%% Plotting 3
%%%%%%%%%%%%%%%%%%%%%%%
% read data
fgmorig = spdfcdfread('/Users/Andrea/Documents/MSci Project/Data/fgm/mms1_fgm_brst_l2_20151127013434_v4.18.0.cdf',... 
    'ShowProgress',true, 'KeepEpochAsIs',true);
% Magnetic field (B) 
B(:,1:4) = fgmorig{1,2}; %Bx,By,Bz,B (GSE- in nT)
% Spacecraft Position
r(:,1:4) = fgmorig{1,8}; %x, y, z, r (GSE- in m)

%%%%%%%%%%%%%%%%%%%%%%%
% specify time interval
tint=[irf_time([2016 01 21 01 06 41]) irf_time([2016 01 21 01 06 52])]; % time interval

%%%%%%%%%%%%%%%%%%%%%%%
% initialize figure
h=irf_plot(3,'newfigure'); % 3 subplots, remove middle (fast solution)
delete(h(2));				% remove middle panel
h(2)=[];					% remove handle
set(h(2),'position',get(h(2),'position')+[0.2 0.1 -0.2 0]); % move panel 2 up a bit

%%%%%%%%%%%%%%%%%%%%%%%
% top panel
hca=h(1);
% read data
%B=irf_get_data('B_vec_xyz_gse__C1_CP_FGM_5VPS','caa','mat');
tirf=irf_time(fgmorig{1,1}, 'ttns>epoch'); % time
tstate=irf_time(fgmorig{1,7}, 'ttns>epoch'); % time for the spacecraft positions
Bwtime(:,1) = tirf;
Bwtime(:,2:4) = B(:,1:3);
% plot
irf_plot(hca,Bwtime);
ylabel(hca,'B [nT] GSE');
%irf_zoom(hca,'y',[-50 50]);
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05]);
irf_legend(hca,{'MMS3'},[0.98 0.98],'color','k');

%%%%%%%%%%%%%%%%%%%%%%%
% top panel
hca=h(2);
% plot
irf_plot(hca,Bwtime);
ylabel(hca,'B [nT] GSE');
%irf_zoom(hca,'y',[-50 50])
irf_legend(hca,{'B_X','B_Y','B_Z'},[0.98 0.05])
irf_legend(hca,{'MMS3'},[0.98 0.98],'color','k')

%%%%%%%%%%%%%%%%%%%%%%%%
% changes to all figure
irf_pl_number_subplots(h);
%irf_zoom(h(1),'x',tint);
% zoom figure 2 to smaller interval
tzoom=tint(1)+diff(tint)*[0.6 0.8];
irf_zoom(h(2),'x',tint);
R1(:,1) = tstate;
R1(:,2:4)= r(:,1:3);
Units= irf_units; % to get RE value
R1RE=irf_tappl(R1,'/Units.RE*Units.km'); % 
xx=get(gcf,'userdata');
tst=xx.t_start_epoch;
xlab={'X (RE)','Y (RE)','Z (RE)'};
irf_timeaxis(h(1),tst,R1RE,xlab);
irf_timeaxis(h(1),'nodate');
% connect zoom in region
irf_plot_zoomin_lines_between_panels(h(1),h(2));
%
%irf_legend(h(1),'Example 3',[1.0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%
% add interval mark
irf_pl_mark(h(1),tint,'yellow')
irf_pl_mark(h(2),tzoom,'blue')