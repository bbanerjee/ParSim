ZRPoint = [2.98457 8.0674 0];
ZRTable = [-3.97635e+06 -257.694 0; -1.98529e+06 -182.085 0; 5773.5 0 0; -1.98529e+06 182.085 0; -3.97635e+06 257.694 0; -5.96741e+06 315.686 0; -7.95847e+06 364.566 0; -9.94954e+06 407.627 0; -1.19406e+07 446.555 0; -1.39317e+07 482.351 0; -1.59227e+07 515.668 0; -1.79138e+07 546.96 0; -1.99048e+07 576.556 0; -2.18959e+07 604.705 0; -2.3887e+07 631.6 0; -2.5878e+07 657.397 0; -2.78691e+07 682.218 0; -2.98602e+07 706.168 0; -3.18512e+07 729.331 0; -3.38423e+07 751.781 0; -3.58333e+07 773.58 0; -3.78244e+07 794.782 0; -3.98155e+07 815.432 0; -4.18065e+07 835.572 0; -4.37976e+07 855.238 0; -4.57887e+07 874.462 0; -4.77797e+07 893.272 0; -4.97708e+07 911.694 0; -5.17618e+07 929.751 0; -5.37529e+07 947.464 0; -5.5744e+07 964.852 0; -5.7735e+07 981.932 0; -5.79341e+07 983.64 0; -5.7954e+07 983.811 0];
ZRSpline = [-989758 -91.0424 0; -925601 -84.9729 0; -865869 -78.9034 0; -810562 -72.8339 0; -759679 -66.7644 0; -713221 -60.6949 0; -671188 -54.6254 0; -633579 -48.5559 0; -600394 -42.4864 0; -571635 -36.417 0; -547299 -30.3475 0; -527389 -24.278 0; -511903 -18.2085 0; -500841 -12.139 0; -494204 -6.06949 0; -491992 0 0; -494204 6.06949 0; -500841 12.139 0; -511903 18.2085 0; -527389 24.278 0; -547299 30.3475 0; -571635 36.417 0; -600394 42.4864 0; -633579 48.5559 0; -671188 54.6254 0; -713221 60.6949 0; -759679 66.7644 0; -810562 72.8339 0; -865869 78.9034 0; -925601 84.9729 0; -989758 91.0424 0];
ZRClose = [-491992 0 0];
%plot(ZRPoint(1), ZRPoint(2), 'kx'); hold on;
%plot(ZRTable(:,1), ZRTable(:,2), 'b-'); hold on;
%plot(ZRSpline(:,1), ZRSpline(:,2), 'r-'); hold on;
%plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'g-'); hold on;

ZRPoint = [5196.15 0 0]
ZRTable = [-17.3205 -173.205 0; 16.9741 -1.73205 0; 17.3205 0 0; 16.9741 1.73205 0; -17.3205 173.205 0; -692.82 866.025 0; -1385.64 1039.23 0; -2771.28 1212.44 0; -5542.56 1385.64 0; -11085.1 1558.85 0; -11639.4 1576.17 0; -11694.8 1577.9 0] 
ZRSpline = [17.1473 -0.866025 0; 17.1585 -0.80829 0; 17.1689 -0.750555 0; 17.1785 -0.69282 0; 17.1873 -0.635085 0; 17.1954 -0.57735 0; 17.2027 -0.519615 0; 17.2093 -0.46188 0; 17.215 -0.404145 0; 17.22 -0.34641 0; 17.2243 -0.288675 0; 17.2277 -0.23094 0; 17.2304 -0.173205 0; 17.2324 -0.11547 0; 17.2335 -0.057735 0; 17.2339 0 0; 17.2335 0.057735 0; 17.2324 0.11547 0; 17.2304 0.173205 0; 17.2277 0.23094 0; 17.2243 0.288675 0; 17.22 0.34641 0; 17.215 0.404145 0; 17.2093 0.46188 0; 17.2027 0.519615 0; 17.1954 0.57735 0; 17.1873 0.635085 0; 17.1785 0.69282 0; 17.1689 0.750555 0; 17.1585 0.80829 0; 17.1473 0.866025 0] 
ZRClose = [17.2339 0 0]

%plot(ZRPoint(1), ZRPoint(2), 'kx'); hold on;
%plot(ZRTable(:,1), ZRTable(:,2), 'b-'); hold on;
%plot(ZRSpline(:,1), ZRSpline(:,2), 'r-'); hold on;
%plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'g-'); hold on;

ZRPoint = [5196.15 1732.05 0]
ZRTable = [-17.3205 -173.205 0; 16.9741 -1.73205 0; 17.3205 0 0; 16.9741 1.73205 0; -17.3205 173.205 0; -692.82 866.025 0; -1385.64 1039.23 0; -2771.28 1212.44 0; -5542.56 1385.64 0; -11085.1 1558.85 0; -11639.4 1576.17 0; -11694.8 1577.9 0] 
ZRSpline = [-0.173205 87.4686 0; -1.67258 93.474 0; -3.88441 100.059 0; -6.80869 107.223 0; -10.4454 114.966 0; -14.7946 123.288 0; -19.8562 132.19 0; -25.6303 141.671 0; -32.1168 151.731 0; -39.3158 162.371 0; -47.2273 173.59 0; -55.8511 185.388 0; -65.1875 197.766 0; -75.2362 210.722 0; -85.9975 224.258 0; -97.4712 238.373 0; -109.657 253.068 0; -122.556 268.342 0; -136.167 284.195 0; -150.49 300.627 0; -165.526 317.639 0; -181.275 335.23 0; -197.736 353.4 0; -214.909 372.149 0; -232.795 391.478 0; -251.393 411.386 0; -270.703 431.873 0; -290.726 452.94 0; -311.462 474.586 0; -332.91 496.811 0; -355.07 519.615 0] 
ZRClose = [-1.67258 93.474 0]

plot(ZRPoint(1), ZRPoint(2), 'kx'); hold on;
plot(ZRTable(:,1), ZRTable(:,2), 'b-'); hold on;
plot(ZRSpline(:,1), ZRSpline(:,2), 'r-'); hold on;
plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'g-'); hold on;