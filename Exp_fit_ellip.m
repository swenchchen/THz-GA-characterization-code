function output=Exp_fit_ellip(read_file,write_file,phase_in)
% clear all

% load data
load(read_file);
saved_file=write_file;
d_sc=20e-6;  % stratum corneum thickness

% phase clibration searching range
phase_Bre_range=[phase_in.phase_Bre-8,phase_in.phase_Bre+8]; 
phase_ATR_range=[phase_in.phase_ATR-8,phase_in.phase_ATR+8];

%% pre setting
% cut fitting range
start_ind=1;
stop_length=find(freq_search<=1);
stop_ind=numel(stop_length);

freq_search=freq_search(start_ind:stop_ind);
ATR_p_ratio_exp=ATR_p_ratio_exp(start_ind:stop_ind);
ATR_s_ratio_exp=ATR_s_ratio_exp(start_ind:stop_ind);
Bre_p_ratio_exp=Bre_p_ratio_exp(start_ind:stop_ind);
Bre_s_ratio_exp=Bre_s_ratio_exp(start_ind:stop_ind);
water_nk=water_dd_vib_model(freq_search,30);

% constant parameters
c=299792458;
n_si=3.418-1i*0.0001;
n_air=1-1i*0.0001;
theta_ATR=56.94; % i.e. thetai1 in the paper
theta_Bre=33.06; % i.e. thetai2 in the paper

% GA parameters
num_pop=10000;  % number of population, must be even
num_itp=50; % number of iteration

n_sc_p=[1.3,2.25];
k_sc_p=[0.1,0.65];  % nk of the sc-e component 
n_sc_s=[1.35,2.25];
k_sc_s=[0.15,1]  ;% nk of the sc-o component 
n_ep=[1.8,3];
k_ep=[0.3,1.4];  % nk of the ep component 

phase_ATR_shift_range=(phase_ATR_range)*pi/180;
phase_Bre_shift_range=(phase_Bre_range)*pi/180;

%%
% calculate the theoretical reference 
Bre_p_ref_theory=r_trilayer(n_si,n_air,n_air,0,theta_Bre,freq_search,'p');
Bre_s_ref_theory=r_trilayer(n_si,n_air,n_air,0,theta_Bre,freq_search,'s');
ATR_p_ref_theory=r_trilayer(n_si,n_air,n_air,0,theta_ATR,freq_search,'p');
ATR_s_ref_theory=r_trilayer(n_si,n_air,n_air,0,theta_ATR,freq_search,'s');

%% initialization
for r=1:3 % if multiple runs are desired, r>1
    nk_sc_p_gene=zeros(num_pop,numel(freq_search));
    nk_sc_s_gene=zeros(num_pop,numel(freq_search));
    nk_ep_gene=zeros(num_pop,numel(freq_search));
    % assign nk by exponential functions
    parfor p=1:num_pop
        nk_sc_p_gene(p,:)=exp_assign(freq_search,n_sc_p,k_sc_p,water_nk,60);
        nk_sc_s_gene(p,:)=exp_assign(freq_search,n_sc_s,k_sc_s,water_nk,70);
        nk_ep_gene(p,:)=exp_assign(freq_search,n_ep,k_ep,water_nk,100);
    end
    phase_ATR_gene=rand(1,num_pop)*(phase_ATR_shift_range(2)-phase_ATR_shift_range(1))+phase_ATR_shift_range(1);
    phase_Bre_gene=rand(1,num_pop)*(phase_Bre_shift_range(2)-phase_Bre_shift_range(1))+phase_Bre_shift_range(1);

    % GA starts
    error_itp=linspace(0,0,num_itp);
    phase_ATR_gene_best=linspace(0,0,num_itp);
    phase_Bre_gene_best=linspace(0,0,num_itp);
    % iteration start
    for x=1:num_itp
        parfor p=1:num_pop
            % the corresponding sample reflections
            ATR_p_theory=r_trilayer_uniaxial(n_si,nk_sc_s_gene(p,:),nk_sc_p_gene(p,:),nk_ep_gene(p,:),d_sc,theta_ATR,freq_search,'p');
            ATR_s_theory=r_trilayer_uniaxial(n_si,nk_sc_s_gene(p,:),nk_sc_p_gene(p,:),nk_ep_gene(p,:),d_sc,theta_ATR,freq_search,'s');
            Bre_p_theory=r_trilayer_uniaxial(n_si,nk_sc_s_gene(p,:),nk_sc_p_gene(p,:),nk_ep_gene(p,:),d_sc,theta_Bre,freq_search,'p');
            Bre_s_theory=r_trilayer_uniaxial(n_si,nk_sc_s_gene(p,:),nk_sc_p_gene(p,:),nk_ep_gene(p,:),d_sc,theta_Bre,freq_search,'s');
            % the corresponding reflection ratios by comparing to the references
            ATR_p_ratio_theory(p,:)=ATR_p_theory./ATR_p_ref_theory;
            ATR_s_ratio_theory(p,:)=ATR_s_theory./ATR_s_ref_theory;
            Bre_p_ratio_theory(p,:)=Bre_p_theory./Bre_p_ref_theory;
            Bre_s_ratio_theory(p,:)=Bre_s_theory./Bre_s_ref_theory;
            % error compared to the experimental data
            ATR_p_abs_error=1*sum((abs(ATR_p_ratio_theory(p,:)-ATR_p_ratio_exp.*exp(1i*freq_search*(phase_ATR_gene(p))))./abs(ATR_p_ratio_exp)).^2);
            ATR_s_abs_error=1*sum((abs(ATR_s_ratio_theory(p,:)-ATR_s_ratio_exp.*exp(1i*freq_search*(phase_ATR_gene(p))))./abs(ATR_s_ratio_exp)).^2);
            Bre_p_abs_error=1*sum((abs(Bre_p_ratio_theory(p,:)-Bre_p_ratio_exp.*exp(1i*freq_search*(phase_Bre_gene(p))))./abs(Bre_p_ratio_exp)).^2);
            Bre_s_abs_error=1*sum((abs(Bre_s_ratio_theory(p,:)-Bre_s_ratio_exp.*exp(1i*freq_search*(phase_Bre_gene(p))))./abs(Bre_s_ratio_exp)).^2);
            % evaluation function
            error_pop(p)=sqrt(mean(ATR_p_abs_error+ATR_s_abs_error+Bre_p_abs_error+Bre_s_abs_error)); % newly added to weight
        end
         % rank from low to high
        [error_pop_sort,pop_sort_ind]=sort(error_pop);
        error_itp(x)=error_pop_sort(1);

        nk_sc_p_pop_sort=nk_sc_p_gene(pop_sort_ind(1:end/2+1),:);
        nk_sc_s_pop_sort=nk_sc_s_gene(pop_sort_ind(1:end/2+1),:);
        nk_ep_pop_sort=nk_ep_gene(pop_sort_ind(1:end/2+1),:);
        phase_ATR_gene_sort=phase_ATR_gene(pop_sort_ind(1:end/2+1));
        phase_Bre_gene_sort=phase_Bre_gene(pop_sort_ind(1:end/2+1));
        % record the best results achieved so far
        ATR_p_ratio_theory_best=ATR_p_ratio_theory(pop_sort_ind(1),:);
        ATR_s_ratio_theory_best=ATR_s_ratio_theory(pop_sort_ind(1),:);
        Bre_p_ratio_theory_best=Bre_p_ratio_theory(pop_sort_ind(1),:);
        Bre_s_ratio_theory_best=Bre_s_ratio_theory(pop_sort_ind(1),:);

        nk_sc_p_best=nk_sc_p_pop_sort(1,:);
        nk_sc_s_best=nk_sc_s_pop_sort(1,:);
        nk_ep_best=nk_ep_pop_sort(1,:);
        phase_ATR_gene_best(x)=phase_ATR_gene_sort(1);
        phase_Bre_gene_best(x)=phase_Bre_gene_sort(1);

        % crossover
        nk_sc_p_new_born=(nk_sc_p_pop_sort(1:end-1,:)+nk_sc_p_pop_sort(2:end,:))/2;
        nk_sc_s_new_born=(nk_sc_s_pop_sort(1:end-1,:)+nk_sc_s_pop_sort(2:end,:))/2;
        nk_ep_new_born=(nk_ep_pop_sort(1:end-1,:)+nk_ep_pop_sort(2:end,:))/2;
        phase_ATR_gene_new_born=(phase_ATR_gene_sort(1:end-1)+phase_ATR_gene_sort(2:end))/2;
        phase_Bre_gene_new_born=(phase_Bre_gene_sort(1:end-1)+phase_Bre_gene_sort(2:end))/2;

        nk_sc_p_gene=[nk_sc_p_pop_sort(1:end-1,:);nk_sc_p_new_born];
        nk_sc_s_gene=[nk_sc_s_pop_sort(1:end-1,:);nk_sc_s_new_born];
        nk_ep_gene=[nk_ep_pop_sort(1:end-1,:);nk_ep_new_born];
        phase_ATR_gene=[phase_ATR_gene_sort(1:end-1),phase_ATR_gene_new_born];
        phase_Bre_gene=[phase_Bre_gene_sort(1:end-1),phase_Bre_gene_new_born];

        % mutation
        if x<=num_itp/4
            mut_rate=0.3;
            mut_step=0.15;
        elseif x>num_itp/4 && x<=num_itp*2/4
            mut_rate=0.4;
            mut_step=0.1;
        elseif x>num_itp*2/4 && x<=num_itp*3/4
            mut_rate=0.55;
            mut_step=0.08;
        else
            mut_rate=0.7;
            mut_step=0.04;
        end
        mut_flag=find(rand(1,num_pop/2)<mut_rate)+num_pop/2;
        if ~isempty(mut_flag)
            nk_sc_p_mut_list=nk_sc_p_gene(mut_flag,:);
            nk_sc_s_mut_list=nk_sc_s_gene(mut_flag,:);
            nk_ep_mut_list=nk_ep_gene(mut_flag,:);
            phase_ATR_mut_list=phase_ATR_gene(mut_flag);
            phase_Bre_mut_list=phase_Bre_gene(mut_flag);
            parfor m=1:numel(mut_flag)
                mut_para_flag=rand*3.6;
                if mut_para_flag<0.7 % mut nk_sc_p
                    nk_sc_p_mut_list(m,:)=exp_mutate(freq_search,nk_sc_p_mut_list(m,:),mut_step);
                elseif mut_para_flag<2.5 % mut nk_sc_s
                    nk_sc_s_mut_list(m,:)=exp_mutate(freq_search,nk_sc_s_mut_list(m,:),mut_step);
                elseif mut_para_flag<3  % mut nk_ep
                    nk_ep_mut_list(m,:)=exp_mutate(freq_search,nk_ep_mut_list(m,:),mut_step);
                elseif mut_para_flag<3.2 % mutate phase shift
                    phase_Bre_mut_list(m)=(1+(rand*2-1)*mut_step)*phase_Bre_mut_list(m);
                else
                    phase_ATR_mut_list(m)=(1+(rand*2-1)*mut_step)*phase_ATR_mut_list(m);
                end
            end
            nk_sc_p_gene(mut_flag,:)=nk_sc_p_mut_list;
            nk_sc_s_gene(mut_flag,:)=nk_sc_s_mut_list;
            nk_ep_gene(mut_flag,:)=nk_ep_mut_list;
            phase_ATR_gene(mut_flag)=phase_ATR_mut_list;
            phase_Bre_gene(mut_flag)=phase_Bre_mut_list;
        end
        % result plot
        if x==1 && r==1
            figure('units','normalized','outerposition',[0 0 1 1])

        end
        subplot(2,3,1)
        plot(1:numel(error_itp),(error_itp))
        hold off
        set(gca,'fontsize',14)
        legend(['loop',int2str(x),' of run',int2str(r)])

        subplot(2,3,4)
        plot(1:numel(phase_ATR_gene_best),phase_ATR_gene_best*180/pi)
        hold all
        plot(1:numel(phase_Bre_gene_best),phase_Bre_gene_best*180/pi)
        hold off
        set(gca,'fontsize',14)

        subplot(2,3,2)
        plot(freq_search,180/pi*angle(ATR_p_ratio_theory_best),'b')
        hold all
        plot(freq_search,180/pi*(angle(ATR_p_ratio_exp)+freq_search.*phase_ATR_gene_best(x)),'b*')

        plot(freq_search,180/pi*angle(ATR_s_ratio_theory_best),'g')
        plot(freq_search,180/pi*(angle(ATR_s_ratio_exp)+freq_search.*phase_ATR_gene_best(x)),'g*')

        plot(freq_search,180/pi*angle(Bre_p_ratio_theory_best),'r')
        plot(freq_search,180/pi*(angle(Bre_p_ratio_exp)+freq_search.*phase_Bre_gene_best(x)),'r*')

        plot(freq_search,180/pi*angle(Bre_s_ratio_theory_best),'k')
        plot(freq_search,180/pi*(angle(Bre_s_ratio_exp)+freq_search.*phase_Bre_gene_best(x)),'k*')
        hold off
        set(gca,'fontsize',14)
        xlim([0.1,1])

        subplot(2,3,5)
        plot(freq_search,abs(ATR_p_ratio_theory_best),'b')
        hold all
        plot(freq_search,abs(ATR_p_ratio_exp),'b*')

        plot(freq_search,abs(ATR_s_ratio_theory_best),'g')
        plot(freq_search,abs(ATR_s_ratio_exp),'g*')

        plot(freq_search,abs(Bre_p_ratio_theory_best),'r')
        plot(freq_search,abs(Bre_p_ratio_exp),'r*')

        plot(freq_search,abs(Bre_s_ratio_theory_best),'k')
        plot(freq_search,abs(Bre_s_ratio_exp),'k*')
        hold off
        set(gca,'fontsize',14)
        xlim([0.1,1])

        subplot(2,3,3)
        plot(freq_search,real(nk_sc_p_best),'b-','linewidth',2)
        hold all
        plot(freq_search,real(nk_sc_s_best),'r-','linewidth',2)
        plot(freq_search,real(nk_ep_best),'k-','linewidth',2)
        plot(freq_search,real(water_nk),'--','color',[0.2,0.2,0.2])
        hold off
        set(gca,'fontsize',14)

        subplot(2,3,6)
        plot(freq_search,-400*pi*freq_search/3.*imag(water_nk),'--','color',[0.2,0.2,0.2])
        hold all
        plot(freq_search,-400*pi*freq_search/3.*imag(nk_sc_p_best),'b-','linewidth',2)
        plot(freq_search,-400*pi*freq_search/3.*imag(nk_sc_s_best),'r-','linewidth',2)
        plot(freq_search,-400*pi*freq_search/3.*imag(nk_ep_best),'k-','linewidth',2)
        hold off
        set(gca,'fontsize',14)
        grid('on')
        pause(0.02)
    end
    % output the optimized results if multiple runs are applied
    error_itp_r(r,:)=error_itp;
    ATR_p_ratio_theory_best_r(r,:)=ATR_p_ratio_theory_best;
    ATR_s_ratio_theory_best_r(r,:)=ATR_s_ratio_theory_best;
    Bre_p_ratio_theory_best_r(r,:)=Bre_p_ratio_theory_best;
    Bre_s_ratio_theory_best_r(r,:)=Bre_s_ratio_theory_best;
    phase_ATR_gene_best_r(r,:)=phase_ATR_gene_best;
    phase_Bre_gene_best_r(r,:)=phase_Bre_gene_best;
    nk_sc_p_best_r(r,:)=nk_sc_p_best;
    nk_sc_s_best_r(r,:)=nk_sc_s_best;
    nk_ep_best_r(r,:)=nk_ep_best;
    
end
%%
% final optimized results by averaging the best results in multiple runs
error_itp=mean(error_itp_r);
ATR_p_ratio_theory_best=mean(ATR_p_ratio_theory_best_r);
ATR_s_ratio_theory_best=mean(ATR_s_ratio_theory_best_r);
Bre_p_ratio_theory_best=mean(Bre_p_ratio_theory_best_r);
Bre_s_ratio_theory_best=mean(Bre_s_ratio_theory_best_r);
phase_Bre_gene_best=mean(phase_Bre_gene_best_r(:,end));
phase_ATR_gene_best=mean(phase_ATR_gene_best_r(:,end));
nk_sc_p_best=mean(nk_sc_p_best_r);
nk_sc_s_best=mean(nk_sc_s_best_r);
nk_ep_best=mean(nk_ep_best_r);

% output
output.phase_Bre=phase_Bre_gene_best*180/pi;
output.phase_ATR=phase_ATR_gene_best*180/pi;
output.nkpsc=nk_sc_p_best;
output.nkssc=nk_sc_s_best;
output.nkep=nk_ep_best;

% save data
save(saved_file,'error_itp','ATR_p_ratio_theory_best','ATR_s_ratio_theory_best','Bre_p_ratio_theory_best','Bre_s_ratio_theory_best','phase_ATR_gene_best','phase_Bre_gene_best','water_nk','nk_sc_p_best','nk_sc_s_best','nk_ep_best','freq_search');