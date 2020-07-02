function output=EMT_LLL(n1,n2,f2_array)

% n1 n2 could be frequency dependent
% f2 is an array related to the layer distribution

if mean(f2_array)>1 % it's in percentage format
    f2_array=f2_array/100;
end

epsilon1=n1.^2;
epsilon2=n2.^2;
for k=1:numel(f2_array)
    f2_c=f2_array(k);
    epsilon_LLL(k,:)=((epsilon2.^(1/3)-epsilon1.^(1/3)).*f2_c+epsilon1.^(1/3)).^3;
end

n_eff=sqrt(epsilon_LLL);

output=n_eff;