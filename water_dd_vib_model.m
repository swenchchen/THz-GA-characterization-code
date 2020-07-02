function water_nk=water_dd_vib_model(freq,T)

a1=80.6972;   b1=0.00436598;    c1=0.13673;   d1=693.48;
a2=3.5997;   b2=0.017544;   c2=0.035735;   d2=356.0099;
tc=139.377;

epsilon_s=87.9144-0.404399.*T+9.58726*10^-4.*T.^2-1.32892*10^-6.*T.^3;

delta1=a1.*exp(-b1.*T);
tao1=c1.*exp(d1./(T+tc));

delta2=a2.*exp(-b2.*T);
tao2=c2.*exp(d2./(T+tc));

relaxation_r1=(2.*pi.*freq).^2.*(tao1.^2.*delta1./(1+(2.*pi.*freq.*tao1).^2));
relaxation_r2=(2.*pi.*freq).^2.*(tao2.^2.*delta2./(1+(2.*pi.*freq.*tao2).^2));

relaxation_i1=tao1.*delta1./(1+(2.*pi.*freq.*tao1).^2);
relaxation_i2=tao2.*delta2./(1+(2.*pi.*freq.*tao2).^2);

p0=0.8379; p1=-0.006118;  p2=-0.000012936;
p3=4.235901;  p4=-0.01426;  p5=0.000273815;   p6=-1.2469e-5;
p7=9.618e-2;  p8=1.79578e-4;  p9=-9.31e-6;  p10=1.655e-7;
delta4=p0+p1.*T+p2.*T.^2;
f0=p3+p4.*T+p5.*T.^2+p6.*T.^3;
tao4=p7+p8.*T+p9.*T.^2+p10.*T.^3;
resonance_r1=-(2*pi*tao4).^2.*delta4/2.*((freq.*(f0+freq))./(1+(2*pi.*tao4.*(f0+freq)).^2) - (freq.*(f0-freq))./(1+(2*pi.*tao4.*(f0-freq)).^2));
resonance_i1=(pi.*freq.*delta4.*tao4).*(1./(1+(2*pi.*tao4.*(f0+freq)).^2) + 1./(1+(2*pi.*tao4.*(f0-freq)).^2));

epsilon_r=epsilon_s-(relaxation_r1+relaxation_r2)+resonance_r1;
epsilon_i=(2.*pi.*freq).*(relaxation_i1+relaxation_i2)+resonance_i1;

epsilon=epsilon_r-1i*epsilon_i;
water_nk=sqrt(epsilon);