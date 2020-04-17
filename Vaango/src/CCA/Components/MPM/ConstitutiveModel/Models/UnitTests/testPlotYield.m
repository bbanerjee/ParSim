line = [[10 0 0];[-10 100 0]];
seg = [[10 0 0];[9.8 1 0];[-10 100 0]];
pt = [2 39 0]; cspline = [2.19231 39.0385 0]; cline = [2.19231 39.0385 0];
%plot(line(:,1), line(:,2), 'b-'); hold on;
%plot(seg(:,1), seg(:,2), 'g-');
%plot(pt(1), pt(2), 'rx');
%plot(cspline(1), cspline(2), 'ko');
%plot(cline(1), cline(2), 'ms');
%axis equal;

line = [[-400 500 0];[-1397 674.625 0]];
seg = [[-400 500 0];[-800 600 0];[-1397 674.625 0]];
pt = [-1000 625 0]; cspline = [-1000.1 624.294 0]; cline = [-1003.38 605.683 0];
%plot(line(:,1), line(:,2), 'b-'); hold on;
%plot(seg(:,1), seg(:,2), 'g-');
%plot(pt(1), pt(2), 'rx');
%plot(cspline(1), cspline(2), 'ko');
%plot(cline(1), cline(2), 'ms');
axis equal;

line = [[-400 500 0]; [-1397 674.625 0] ];
seg = [[-400 500 0];[-800 600 0];[-1397 674.625 0]];
pt = [-1000 605 0];
cspline = [-997.343 623.907 0]; 
cline = [-999.985 605.088 0];
%plot(line(:,1), line(:,2), 'b-'); hold on;
%plot(seg(:,1), seg(:,2), 'g-');
%plot(pt(1), pt(2), 'rx');
%plot(cspline(1), cspline(2), 'ko');
%plot(cline(1), cline(2), 'ms');
%axis equal;

line = [[-1651.84 637.352 0]; [-1742.87 580.721 0] ];
seg = [[-1651.84 637.352 0];[-1698.5 611.549 0];[-1742.87 580.721 0]];
pt = [-1700 610.613 0];
cspline = [-1699.64 610.032 0]; 
cline = [-1698.55 608.289 0];
%plot(line(:,1), line(:,2), 'b-'); hold on;
%plot(seg(:,1), seg(:,2), 'g-');
%plot(pt(1), pt(2), 'rx');
%plot(cspline(1), cspline(2), 'ko');
%plot(cline(1), cline(2), 'ms');
%axis equal;

line = [[-9597.67 499.964 0];[-9818.9 344.36 0] ];
seg = [[-9597.67 499.964 0];[-9718.64 424.185 0];[-9818.9 344.36 0]];
pt = [-9700 437.045 0];
cspline = [-9698.95 435.494 0];
cline = [-9695.74 430.987 0];
%plot(line(:,1), line(:,2), 'b-'); hold on;
%plot(seg(:,1), seg(:,2), 'g-');
%plot(pt(1), pt(2), 'rx');
%plot(cspline(1), cspline(2), 'ko');
%plot(cline(1), cline(2), 'ms');
%axis equal;

ZRSeg0 = [42435.24478543749 0 0];
ZRSeg1 = [41539.19694987215 955.9687004835084 0];
ZRSeg2 = [-37385.81438664001 85158.88458098272 0];
ZRPoint = [11497.90431914181 33223.01418229147 0];
ZRClose = [11385.95668977881 33125.62835779594 0];
%plot([ZRSeg0(1) ZRSeg1(1) ZRSeg2(1)],[ZRSeg0(2) ZRSeg1(2) ZRSeg2(2)],'r-','linewidth',2); hold on;
plot(ZRPoint(1), ZRPoint(2),'bx','linewidth',2); hold on;
plot(ZRClose(1), ZRClose(2),'mx','linewidth',2); hold on;

ZRSeg0 = [42435.24478543749 0 0];
ZRSeg1 = [41539.19694987215 955.9687004835084 0];
ZRSeg2 = [-37385.81438664001 85158.88458098272 0];
ZRPoint = [11489.70691392369 33231.75106159761 0];
ZRClose = [9671.033650183563 34955.23240851802 0];

plot([ZRSeg0(1) ZRSeg1(1) ZRSeg2(1)],[ZRSeg0(2) ZRSeg1(2) ZRSeg2(2)],'gx--','linewidth',2); hold on;
plot(ZRPoint(1), ZRPoint(2),'o','color', [0.25 0.75 0.15],'linewidth',2); hold on;
plot(ZRClose(1), ZRClose(2),'o','color', [0.75 0.25 0.15],'linewidth',2); hold on;

