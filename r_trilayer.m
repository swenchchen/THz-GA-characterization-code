function output=r_trilayer(n1,n2,n3,thickness,theta_i,freq,polarization)
%%
c=299792458;
theta_i=theta_i*pi/180;
theta_t12=asin(n1.*sin(theta_i)./n2);
theta_t23=asin(n1.*sin(theta_i)./n3);

if polarization=='s'
    r12=(n1.*cos(theta_i)-n2.*cos(theta_t12))./(n1.*cos(theta_i)+n2.*cos(theta_t12));
    r23=(n2.*cos(theta_t12)-n3.*cos(theta_t23))./(n2.*cos(theta_t12)+n3.*cos(theta_t23));
elseif polarization=='p'
    r12=(n2.*cos(theta_i)-n1.*cos(theta_t12))./(n2.*cos(theta_i)+n1.*cos(theta_t12));
    r23=(n3.*cos(theta_t12)-n2.*cos(theta_t23))./(n3.*cos(theta_t12)+n2.*cos(theta_t23));
else
    msg='polarization input error'
end

beta=2*pi.*freq*10^12.*thickness.*n2.*cos(theta_t12)./c;

r123=(r12+r23.*exp(-2i*beta))./(1+r12.*r23.*exp(-2i*beta));

output=r123;