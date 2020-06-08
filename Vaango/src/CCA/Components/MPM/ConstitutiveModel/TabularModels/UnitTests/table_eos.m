format long e

v = [0.12 150 0.25]
sal = [0.1 0.2];
s = (v(1) - sal(1))/(sal(2) - sal(1))
%
temp1 = [100 200 300];
t1 = (v(2) - temp1(1))/(temp1(2) - temp1(1))
vol1 = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
press1 = [10 20 30 40 50 60 70 80];
density1 = [1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80];
u1 = (v(3) - vol1(2))/(vol1(3) - vol1(2))
p1 = (1 - u1)*press1(2) + u1*press1(3)
d1 = (1 - u1)*density1(2) + u1*density1(3)
vol2 = [0.15 0.25 0.35 0.45 0.55];
press2 = [100 200 300 400 500];
density2 = [2.100 2.200 2.300 2.400 2.500];
u2 = (v(3) - vol2(1))/(vol2(2) - vol2(1))
p2 = (1 - u2)*press2(1) + u2*press2(2)
d2 = (1 - u2)*density2(1) + u2*density2(2)
pp1 = (1 - t1)*p1 + t1*p2
dd1 = (1 - t1)*d1 + t1*d2
%
temp2 = [0 400];
t2 = (v(2) - temp2(1))/(temp2(2) - temp2(1))
vol3 = [0.1 0.2 0.3 0.8];
press3 = [15 25 35 45];
density3 = [4.15 4.25 4.35 4.45];
u3 = (v(3) - vol3(2))/(vol3(3) - vol3(2))
p3 = (1 - u3)*press3(2) + u3*press3(3)
d3 = (1 - u3)*density3(2) + u3*density3(3)
vol4 = [0.1 0.45 0.65];
press4 = [150 250 350];
density4 = [5.150 5.250 5.350];
u4 = (v(3) - vol4(1))/(vol4(2) - vol4(1))
p4 = (1 - u4)*press4(1) + u4*press4(2)
d4 = (1 - u4)*density4(1) + u4*density4(2)
pp2 = (1 - t2)*p3 + t2*p4
dd2 = (1 - t2)*d3 + t2*d4
%
ppp = (1 - s)*pp1 + s*pp2
ddd = (1 - s)*dd1 + s*dd2

v = [0.12 230 0.5]
sal = [0.1 0.2];
s = (v(1) - sal(1))/(sal(2) - sal(1))
%
temp1 = [100 200 300];
t1 = (v(2) - temp1(2))/(temp1(3) - temp1(2))
vol1 = [0.15 0.25 0.35 0.45 0.55];
press1 = [100 200 300 400 500];
density1 = [2.100 2.200 2.300 2.400 2.500];
u1 = (v(3) - vol1(4))/(vol1(5) - vol1(4))
p1 = (1 - u1)*press1(4) + u1*press1(5)
d1 = (1 - u1)*density1(4) + u1*density1(5)
vol2 = [0.05 0.75];
press2 = [1000 2000];
density2 = [3.1000 3.2000];
u2 = (v(3) - vol2(1))/(vol2(2) - vol2(1))
p2 = (1 - u2)*press2(1) + u2*press2(2)
d2 = (1 - u2)*density2(1) + u2*density2(2)
pp1 = (1 - t1)*p1 + t1*p2
dd1 = (1 - t1)*d1 + t1*d2
%
temp2 = [0 400];
t2 = (v(2) - temp2(1))/(temp2(2) - temp2(1))
vol3 = [0.1 0.2 0.3 0.8];
press3 = [15 25 35 45];
density3 = [4.15 4.25 4.35 4.45];
u3 = (v(3) - vol3(3))/(vol3(4) - vol3(3))
p3 = (1 - u3)*press3(3) + u3*press3(4)
d3 = (1 - u3)*density3(3) + u3*density3(4)
vol4 = [0.1 0.45 0.65];
press4 = [150 250 350];
density4 = [5.150 5.250 5.350];
u4 = (v(3) - vol4(2))/(vol4(3) - vol4(2))
p4 = (1 - u4)*press4(2) + u4*press4(3)
d4 = (1 - u4)*density4(2) + u4*density4(3)
pp2 = (1 - t2)*p3 + t2*p4
dd2 = (1 - t2)*d3 + t2*d4
%
ppp = (1 - s)*pp1 + s*pp2
ddd = (1 - s)*dd1 + s*dd2
