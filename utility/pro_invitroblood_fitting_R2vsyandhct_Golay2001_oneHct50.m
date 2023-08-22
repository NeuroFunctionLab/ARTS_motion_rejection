
x=zeros(24,2);
x(:,1)=[
1
0.805
0.467
0.591
1 
0.745
0.793
0.808
0.573
0.583
0.384
0.484
1
0.77
0.66
0.469
1
0.867
0.69
0.521
0.809
0.585
0.4
1
];

x(:,2)=[
0.459
0.464
0.456
0.464
0.362
0.363
0.36
0.361
0.362
0.363
0.367
0.362
0.556
0.561
0.561
0.574
0.398
0.405
0.409
0.413
0.501
0.507
0.513
0.517
];
y=[
6.8989
8.6663
25.7813
16.7878
5.8203
9.2866
7.9022
7.812
17.069
16.7678
26.4307
21.3272
7.8463
11.1579
14.9416
24.491
6.553
7.5725
12.8434
20.7378
10.5013
19.7919
31.709
8.0625
];
yy=y;
[temp1,resid, jacob]=nlinfit(x, y,'model_R2_vs_y_hct_Golay2001',[10,10,10,10,10,10]);
temp1;

count1=0;
simr2=zeros(71,24);
for y=0.3:0.01:1
    count1=count1+1;
    count2=0;
    for hct=0.35:0.01:0.58
        count2=count2+1;
simr2(count1,count2)=model_R2_vs_y_hct_Golay2001 (temp1, [y,hct]);
    end
end
hct=Hct;
count=0;
for y=0.3:0.001:1
    count=count+1;
    r2_y(count)=model_R2_vs_y_hct_Golay2001 (temp1, [y,hct]);
end
y=0.3:0.001:1;
measured_r2 = 1000./measured_t2;
abs_dif=abs(r2_y-measured_r2);
[temp2, temp3]=min(abs_dif);
measured_oxy=y(temp3);
