function Yv=prox_invitroblood_batch_R2vsyandhct_Golay2001_oneHct50_1(t2,hct)

t2matrix=t2;	
Hct=hct;
t=length(t2matrix);
%t2matrix(1,1)=57.7202;t2matrix(2,1)=58.2798;t2matrix(3,1)=61.9915; t2matrix(4,1)=61.8519;
myoxy=zeros(t,1);
for i=1:t
    for j=1:1
        measured_t2=t2matrix(i,j);
        pro_invitroblood_fitting_R2vsyandhct_Golay2001_oneHct50;
        myoxy(i,j)=measured_oxy*100;
    end
end
Yv=myoxy;   