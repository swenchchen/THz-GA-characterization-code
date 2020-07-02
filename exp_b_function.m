function output=exp_b_function(b,x_array,Y)

output=(exp(b*x_array(3))-exp(b*x_array(1)))./(exp(b*x_array(2))-exp(b*x_array(1)))-Y;