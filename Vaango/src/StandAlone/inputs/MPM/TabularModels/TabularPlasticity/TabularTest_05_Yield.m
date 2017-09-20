minI1 = -1.0e3;
maxI1 =  2.0e3;
I1 = linspace(minI1, maxI1, 30);
J2 = 0.01*maxI1*(I1 - minI1);
SqrtJ2 = sqrt(J2);
plot(I1/3, SqrtJ2, 'k-', 'LineWidth', 2); hold on;
axis equal
grid on

printf('%f, ', I1/3)
printf('\n')
printf('%f, ', SqrtJ2)
printf('\n')
