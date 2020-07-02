function nk_gene=exp_assign(freq,n_bd,k_bd,water_nk,max_wf)
%%
nk_up_limit=EMT_LLL(1.9-1i*0.3,water_nk,max_wf);  % use an EMT to set the upper limit of the searching range

rand_freq=rand*(freq(end)-freq(1)-0.2)+freq(1)+0.1;
freq_i=[freq(1),rand_freq,freq(end)];

n_rand1=rand*(min(n_bd(2),real(nk_up_limit(1)))-n_bd(1))+n_bd(1);
n_rand3=rand*(min(n_rand1-0.02,real(nk_up_limit(end)))-max(n_bd(1),0.62*(n_rand1-1)+1))+max(n_bd(1),0.62*(n_rand1-1)+1);
n_mid_max=0.95*(freq(end)-rand_freq)./(freq(end)-freq(1)).*(n_rand1-n_rand3)+n_rand3;
n_mid_min=0.1*(freq(end)-rand_freq)./(freq(end)-freq(1)).*(n_rand1-n_rand3)+n_rand3;
n_rand2=rand*(n_mid_max-n_mid_min)+n_mid_min;
n_i=[n_rand1,n_rand2,n_rand3];
n_coef=exp_solution(freq_i,n_i);
n_gene=n_coef(1).*exp(n_coef(2).*freq)+n_coef(3);

k_rand1=rand*(min(k_bd(2),-imag(nk_up_limit(1)))-k_bd(1))+k_bd(1);
k_rand3=rand*(min(k_rand1,-imag(nk_up_limit(end)))-max(k_bd(1),0.36*k_rand1))+max(k_bd(1),0.36*k_rand1);
k_mid_max1=0.95*(freq(end)-rand_freq)./(freq(end)-freq(1)).*(k_rand1-k_rand3)+k_rand3;
k_mid_max2=0.9*k_rand3*(freq_i(3))/(freq_i(2));
k_mid_max=min(k_mid_max1,k_mid_max2);
k_mid_min=0.2*(freq(end)-rand_freq)./(freq(end)-freq(1)).*(k_rand1-k_rand3)+k_rand3;
k_rand2=rand*(k_mid_max-k_mid_min)+k_mid_min;
k_i=[k_rand1,k_rand2,k_rand3];
k_coef=exp_solution(freq_i,k_i);
k_gene=k_coef(1).*exp(k_coef(2).*freq)+k_coef(3);

    
nk_gene=n_gene-1i*k_gene;