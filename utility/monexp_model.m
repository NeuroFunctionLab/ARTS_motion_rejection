function result = monexp_model (beta,x)
%Input:
%	beta   --  parameter of the function, must have 2 element
%	x      --  input vector, must have one column
%Output:
%	result --  corresponding function value, 
%             has the same dimension as x
%Hanzhang Lu
%Date: 07/27/2004



s0 = beta(1);
r2_star=beta(2);

x1 = x(:,1);

result = s0*exp(-x1.*r2_star);
