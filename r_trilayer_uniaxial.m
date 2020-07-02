function output=r_trilayer_uniaxial(n1,n2o,n2e,n3,thickness,theta_i,freq,polarization)
%%
c=299792458;
theta_i=theta_i*pi/180;
theta_t23=asin(n1.*sin(theta_i)./n3);

n_sin_sq=n1.^2.*sin(theta_i).^2;

if polarization=='s'
    r12=(n1.*cos(theta_i)-sqrt(n2o.^2-n_sin_sq)) ./ (n1.*cos(theta_i)+sqrt(n2o.^2-n_sin_sq));
    r23=(sqrt(n2o.^2-n_sin_sq)-n3.*cos(theta_t23)) ./ (sqrt(n2o.^2-n_sin_sq)+n3.*cos(theta_t23));
elseif polarization=='p'
    r12=(n2o.*n2e.*cos(theta_i)-n1.*sqrt(n2e.^2-n_sin_sq)) ./ (n2o.*n2e.*cos(theta_i)+n1.*sqrt(n2e.^2-n_sin_sq));
    r23=(n3.*sqrt(n2e.^2-n_sin_sq)-n2o.*n2e.*cos(theta_t23)) ./ (n3.*sqrt(n2e.^2-n_sin_sq)+n2o.*n2e.*cos(theta_t23));
else
    msg='polarization input error'
end

if polarization=='s'
    beta=(2*pi.*freq*10^12.*thickness./c) .* sqrt(n2o.^2-n_sin_sq);
    r123=(r12+r23.*exp(-2i*beta))./(1+r12.*r23.*exp(-2i*beta));
elseif polarization=='p'
    beta=(2*pi.*freq*10^12.*thickness./c) .* ( (n2o./n2e).*sqrt(n2e.^2-n_sin_sq) );
    r123=(r12+r23.*exp(-2i*beta))./(1+r12.*r23.*exp(-2i*beta));
end

output=r123;