function nk_gene=exp_mutate(freq,nk_gene,mut_step)

section_length=numel(nk_gene)/5;
sampled_ind=round([rand*section_length+1,rand*section_length+2*section_length,rand*section_length+4*section_length-2]);
freq_i=[freq(sampled_ind(1)),freq(sampled_ind(2)),freq(sampled_ind(3))];
n_i=real([nk_gene(sampled_ind(1)),nk_gene(sampled_ind(2)),nk_gene(sampled_ind(3))]);
k_i=-imag([nk_gene(sampled_ind(1)),nk_gene(sampled_ind(2)),nk_gene(sampled_ind(3))]);

mut_flag=round(rand);
if mut_flag==0
    n_coef=exp_solution(freq_i,n_i);
    mut_coef_ind=round(rand*2+1);
    n_coef(mut_coef_ind)=n_coef(mut_coef_ind)*(1+(rand*2-1)*mut_step);
    n_gene=n_coef(1).*exp(n_coef(2).*freq)+n_coef(3);
    k_gene=-imag(nk_gene);
else
    k_coef=exp_solution(freq_i,k_i);
    mut_coef_ind=round(rand*2+1);
    k_coef(mut_coef_ind)=k_coef(mut_coef_ind)*(1+(rand*2-1)*mut_step);
    k_gene=k_coef(1).*exp(k_coef(2).*freq)+k_coef(3);
    n_gene=real(nk_gene);
end

nk_gene=n_gene-1i*k_gene;