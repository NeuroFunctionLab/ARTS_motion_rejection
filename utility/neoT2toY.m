function result = neoT2toY(t2, hct)
% convert T2 values to Yv
% Input: T2: T2 value
%        hct: himatocrit level
% Result: Y: oxygenation value
% Peiying Liu, 06-19-2012

temp1=[-1.0747    1.5148   21.4496   -5.0504   29.4162  242.8936];

 fitted_r2=zeros(1,1001);
 for i=1:1001
     fitted_r2(i)=model_R2_vs_y_hct_Golay2001(temp1, [(i-1)*0.001, hct]);
 end
x=0:0.001:1;

%put the subject's measured t2 value here
measured_t2=t2;

measured_r2=1000/measured_t2;
abs_dif=abs(fitted_r2-measured_r2);
[temp2,temp3]=min(abs_dif);
result=x(temp3);