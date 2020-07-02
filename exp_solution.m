function output_coef=exp_solution(xi,yi)

% x array and y array should be an array including three values
% this function tries to find a,b,c that make a*exp(bx+c)=y fits the input
%%
Y=(yi(3)-yi(1))/(yi(2)-yi(1));
fun=@(b) exp_b_function(b,xi,Y);
b_bd=[-60,-1e-3];
b=fzero(fun,b_bd);

a=(yi(2)-yi(1))./(exp(b*xi(2))-exp(b*xi(1)));

c=yi(1)-a*exp(b*xi(1));

% x_array=xi(1):0.01:xi(3);
% y_array=a.*exp(b.*x_array)+c;

% figure
% plot(xi,yi,'bo')
% hold all
% plot(x_array,y_array,'r','linewidth',2)

output_coef=[a,b,c];
