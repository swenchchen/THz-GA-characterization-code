function output=warmColor(data_num)
pos=[0,25,50,75,100]; % percent position
r_point=[79,246,255,255,242];
g_point=[33,36,25,135,228];
b_point=[211,224,95,46,0];

p_array=0:1:100;
r_array=interp1(pos,r_point,p_array,'linear');
g_array=interp1(pos,g_point,p_array,'linear');
b_array=interp1(pos,b_point,p_array,'linear');

p_used=linspace(0,100,data_num);
r_out=interp1(p_array,r_array,p_used,'linear');
g_out=interp1(p_array,g_array,p_used,'linear');
b_out=interp1(p_array,b_array,p_used,'linear');

output=[r_out',g_out',b_out']/255;