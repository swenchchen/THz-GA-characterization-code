%%
% run GA optimizations for the data measured from 1-31 min
clear all
phase_in.phase_Bre=-4;  phase_in.phase_ATR=-2; % initial phase shift range for the data measured at 1 min
run1=Exp_fit_ellip('Exp_ratio_left_min1.mat','result_left_phase20_min1.mat',phase_in); 
run3=Exp_fit_ellip('Exp_ratio_left_min3.mat','result_left_phase20_min3.mat',run1);
run5=Exp_fit_ellip('Exp_ratio_left_min5.mat','result_left_phase20_min5.mat',run3);
run7=Exp_fit_ellip('Exp_ratio_left_min7.mat','result_left_phase20_min7.mat',run5);
run9=Exp_fit_ellip('Exp_ratio_left_min9.mat','result_left_phase20_min9.mat',run7);
run11=Exp_fit_ellip('Exp_ratio_left_min11.mat','result_left_phase20_min11.mat',run9);
run13=Exp_fit_ellip('Exp_ratio_left_min13.mat','result_left_phase20_min13.mat',run11);
run15=Exp_fit_ellip('Exp_ratio_left_min15.mat','result_left_phase20_min15.mat',run13);
run17=Exp_fit_ellip('Exp_ratio_left_min17.mat','result_left_phase20_min17.mat',run15);
run19=Exp_fit_ellip('Exp_ratio_left_min19.mat','result_left_phase20_min19.mat',run17);
run21=Exp_fit_ellip('Exp_ratio_left_min21.mat','result_left_phase20_min21.mat',run19);
run23=Exp_fit_ellip('Exp_ratio_left_min23.mat','result_left_phase20_min23.mat',run21);
run25=Exp_fit_ellip('Exp_ratio_left_min25.mat','result_left_phase20_min25.mat',run23);
run27=Exp_fit_ellip('Exp_ratio_left_min27.mat','result_left_phase20_min27.mat',run25);
run29=Exp_fit_ellip('Exp_ratio_left_min29.mat','result_left_phase20_min29.mat',run27);
run31=Exp_fit_ellip('Exp_ratio_left_min31.mat','result_left_phase20_min31.mat',run29);

%%
% plot the characterization results
clear all
color_set=warmColor(16);
figure
for k=1:2:31
    load_file=['result_left_phase20_min',int2str(k),'.mat'];
    load(load_file);
    
    phase_ATR_shift((k+1)/2)=phase_ATR_gene_best*180/pi;
    phase_Bre_shift((k+1)/2)=phase_Bre_gene_best*180/pi;
    plot_freq=0.4;
    plot_ind=find(abs(freq_search-plot_freq)==min(abs(freq_search-plot_freq)));
    nk_sc_p((k+1)/2)=nk_sc_p_best(plot_ind);  nk_sc_s((k+1)/2)=nk_sc_s_best(plot_ind);   nk_ep((k+1)/2)=nk_ep_best(plot_ind);
    error_min((k+1)/2)=error_itp(end);
    
    subplot(2,3,1)
    if k==1; plot(freq_search, real(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all
    plot(freq_search,real(nk_sc_p_best),'color',color_set((k+1)/2,:)); 
    xlim([0.1,1]); ylim([1,3.5]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Refractive index'); grid('on');
    subplot(2,3,2)
    if k==1; plot(freq_search, real(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all
    plot(freq_search,real(nk_sc_s_best),'color',color_set((k+1)/2,:)); 
    xlim([0.1,1]); ylim([1,3.5]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Refractive index'); grid('on');
    subplot(2,3,3)
    if k==1; plot(freq_search, real(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all
    plot(freq_search,real(nk_ep_best),'color',color_set((k+1)/2,:)); 
    hold all
    xlim([0.1,1]); ylim([1,3.5]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Refractive index'); grid('on');
    subplot(2,3,4)
    if k==1; plot(freq_search, -400*pi*freq_search/3.*imag(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all;
    plot(freq_search,-400*pi*freq_search/3.*imag(nk_sc_p_best),'color',color_set((k+1)/2,:)); 
    xlim([0.1,1]); ylim([0,250]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Absorption coef. (cm-1)'); grid('on');
    subplot(2,3,5)
    if k==1; plot(freq_search, -400*pi*freq_search/3.*imag(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all;
    plot(freq_search,-400*pi*freq_search/3.*imag(nk_sc_s_best),'color',color_set((k+1)/2,:)); 
    xlim([0.1,1]); ylim([0,250]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Absorption coef. (cm-1)'); grid('on');
    subplot(2,3,6)
    if k==1; plot(freq_search, -400*pi*freq_search/3.*imag(water_nk),'--','color',[0.2,0.2,0.2]); end; hold all;
    plot(freq_search,-400*pi*freq_search/3.*imag(nk_ep_best),'color',color_set((k+1)/2,:));
    hold all
    xlim([0.1,1]); ylim([0,250]); set(gca,'fontsize',18,'fontweight','bold'); xlabel('Frequency (THz)'); ylabel('Absorption coef. (cm-1)'); grid('on');
end

meas_time=1:2:31;
figure
subplot(1,3,1)
plot(meas_time,real(nk_sc_p),'>-','markersize',10,'markerfacecolor','w','linewidth',2)
hold all
plot(meas_time,real(nk_sc_s),'^-','markersize',10,'markerfacecolor','w','linewidth',2)
plot(meas_time,real(nk_ep),'o-','markersize',10,'markerfacecolor','w','linewidth',2)
xlim([0,32]);  set(gca,'fontsize',20,'fontweight','bold'); xlabel('Meas. Time (min)'); ylabel('Refractive index'); grid('on');
subplot(1,3,2)
plot(meas_time,-400*pi*freq_search(plot_ind)/3.*imag(nk_sc_p),'>-','markersize',10,'markerfacecolor','w','linewidth',2)
hold all
plot(meas_time,-400*pi*freq_search(plot_ind)/3.*imag(nk_sc_s),'^-','markersize',10,'markerfacecolor','w','linewidth',2)
plot(meas_time,-400*pi*freq_search(plot_ind)/3.*imag(nk_ep),'o-','markersize',10,'markerfacecolor','w','linewidth',2)
xlim([0,32]);  set(gca,'fontsize',20,'fontweight','bold'); xlabel('Meas. Time (min)'); ylabel('Absorption coef. (cm-1)'); grid('on');
subplot(1,3,3)
plot(meas_time,error_min,'ko-','linewidth',2)
xlim([0,32]);  set(gca,'fontsize',20,'fontweight','bold'); xlabel('Meas. Time (min)'); ylabel('Fitting Error'); grid('on');

figure
for k=1:2:31
    load_file=['Exp_ratio_left_min',int2str(k),'.mat'];
    load(load_file);
    subplot(2,2,1)
    plot(freq_search,180/pi*angle(ATR_p_ratio_exp)+freq_search.*phase_ATR_shift((k+1)/2),'color',color_set((k+1)/2,:))
    set(gca,'fontsize',16);  xlabel('f (THz)');   ylabel('phase (deg.)');  title('ATR-p')
    hold all
    subplot(2,2,2)
    plot(freq_search,180/pi*angle(ATR_s_ratio_exp)+freq_search.*phase_ATR_shift((k+1)/2),'color',color_set((k+1)/2,:))
    set(gca,'fontsize',16);  xlabel('f (THz)');   ylabel('phase (deg.)');  title('ATR-s')
    hold all
    subplot(2,2,3)
    plot(freq_search,180/pi*angle(Bre_p_ratio_exp)+freq_search.*phase_Bre_shift((k+1)/2),'color',color_set((k+1)/2,:))
    set(gca,'fontsize',16);  xlabel('f (THz)');   ylabel('phase (deg.)');  title('Bre-p')
    hold all
    subplot(2,2,4)
    plot(freq_search,180/pi*angle(Bre_s_ratio_exp)+freq_search.*phase_Bre_shift((k+1)/2),'color',color_set((k+1)/2,:))
    set(gca,'fontsize',16);  xlabel('f (THz)');   ylabel('phase (deg.)');  title('Bre-s')
    hold all
end


figure
plot(meas_time,phase_ATR_shift,'o-','linewidth',2)
hold all
plot(meas_time,phase_Bre_shift,'^-','linewidth',2)
set(gca,'fontsize',16);  xlabel('time (min)');   ylabel('phase (deg.)');  legend('ATR phase calibration','Bre phase calibration')