function testElliptic

 phi = pi/4.0;
 mm = 0.2;
 F = EllipticF(phi, mm)

 phi = pi/3.0;
 mm = 0.8;
 F = EllipticF(phi, mm)
 
 phi = pi/2.0;
 mm = 0.8;
 F = EllipticF(phi, mm)
 
%
% Compute Elliptic Integral F(phi|m)
%  0 < phi < pi/2
%  0 < m < 1
%
function [F] = EllipticF(phi, mm)

  % Parameters
  phi_s = 1.249;
  y_s = 0.9;

  % Compute phi_c and m_c
  phi_c = pi/2.0 - phi;
  m_c = 1.0 - mm;

  if (phi < phi_s) 
    F = JacobiSnInv(sin(phi), mm);
  else
    cc = sin(phi_c);
    xx = cc*cc;
    d2 = m_c + mm*xx;
    if (xx < y_s*d2)
      [K, E] = ellipke(m_c);
      F = K - JacobiSnInv(cc/sqrt(d2), mm);
    else
      vv = m_c*(1.0-xx);
      if (vv < x*d2)
        F = JacobiCn(cc, m_c);
      else
        [K, E] = ellipke(m_c);
        F = K - JacobiCnInv(sqrt(vv/d2), m_c);
      end
    end
  end

  return;

%
%  Compute the inverse of the Jacobi sine amplitude function sn(s|m)
%   0 < s < 1
%   0 < m < 1
function [sn_inv] = JacobiSnInv(ss, mm)

  y_a = 0.04095 - 0.00652*mm;
  yy = ss*ss;

  if (yy < y_a)
    sn_inv = ss*SeriesF(yy, mm);
    return;
  end

  pp = 1.0;
  for jj=1:10
    yy = yy/((1.0+sqrt(1.0-yy))*(1.0 + sqrt(1.0-mm*yy)));
    pp = pp*2.0;

    if (yy < y_a)
      sn_inv = pp*sqrt(yy)*SeriesF(yy, mm);
      return;
    end
  end

  err = 'JacobiSnInv: Too many half argument transformations of sn'
  pause

%
%  Compute the inverse of the Jacobi cosine amplitude function cn(s|m)
%   0 < s < 1
%   0 < m < 1
function [cn_inv] = JacobiCnInv(cc, m_c)

  mm = 1.0 - m_c;
  pp = 1.0;
  
  xx = cc*cc;
  for jj=1:10
    if (xx > 0.5)
      cn_inv = pp*JacobiSnInv(sqrt(1.0-xx), mm);
      return;
    end

    dd = sqrt(m_c + mm*xx);
    xx = (sqrt(xx)+dd)/(1.0+dd);
    pp = pp*2.0;
  end

  err = 'JacobiCnInv: Too many half argument transformations of sn'
  pause

%
%  Compute the series expansion used in sine amplitude function sn(s|m)
%   0 < s < 1
%   0 < m < 1
function [series_sum] = SeriesF(yy, mm)

  % Parameters
  u0 = [1.0];
  u1 = [1.0/6.0];
  u2 = [3.0 2.0]/40.0;
  u3 = [5.0 3.0]/112.0;
  u4 = [35.0 20.0 18.0]/1152.0;
  u5 = [63.0 35.0 30.0]/2816.0;
  u6 = [231.0 126.0 105.0 100.0]/13312.0;
  u7 = [429.0 231.0 189.0 175.0]/30720.0;
  u8 = [6435.0 3432.0 2772.0 2520.0 2450.0]/557056.0;
  u9 = [12155.0 6435.0 5148.0 4620.0 4410.0]/1245184.0;
  u10a = [46189.0 24310.0 19305.0 17160.0 16170.0 15876.0]/5505024.0;
  u11 = [88179.0 46189.0 36465.0 32175.0 30030.0 29106.0]/12058624.0;
  u12 = [676039.0 352716.0 277134.0 243100.0 225225.0 216216.0 213444.0]/104857600.0;
  u13 = [1300075.0 676039.0 529074.0 461890.0 425425.0 405405.0 396396.0]/226492416.0;
  u00 = u0(1);
  u10 = u1(1);
  u20 = u2(1); u21 = u2(2);
  u30 = u3(1); u31 = u3(2);
  u40 = u4(1); u41 = u4(2); u42 = u4(3);
  u50 = u5(1); u51 = u5(2); u52 = u5(3);
  u60 = u6(1); u61 = u6(2); u62 = u6(3); u63 = u6(4);
  u70 = u7(1); u71 = u7(2); u72 = u7(3); u73 = u7(4);
  u80 = u8(1); u81 = u8(2); u82 = u8(3); u83 = u8(4); u84 = u8(5);
  u90 = u9(1); u91 = u9(2); u92 = u9(3); u93 = u9(4); u94 = u9(5);
  u100 = u10a(1); u101 = u10a(2); u102 = u10a(3); u103 = u10a(4); u104 = u10a(5); u105 = u10a(6);
  u110 = u11(1); u111 = u11(2); u112 = u11(3); u113 = u11(4); u114 = u11(5); u115 = u11(6);
  u120 = u12(1); u121 = u12(2); u122 = u12(3); u123 = u12(4); u124 = u12(5); u125 = u12(6); u126 = u12(7);
  u130 = u13(1); u131 = u13(2); u132 = u13(3); u133 = u13(4); u134 = u13(5); u135 = u13(6); u136 = u13(7);

  uu0 = u00;
  uu1 = u10 + mm*u10;
  uu2 = u20 + mm*(u21 + mm*u20);
  uu3 = u30 + mm*(u31 + mm*u30);
  uu4 = u40 + mm*(u41 + mm*(u42 + mm*(u42 + (u41 + mm*u40))));
  uu5 = u50 + mm*(u51 + mm*(u52 + mm*(u52 + (u51 + mm*u50))));
  uu6 = u60 + mm*(u61 + mm*(u62 + mm*(u63 + mm*(u63 + mm*(u62 + (u61 + mm*u60))))));
  uu7 = u70 + mm*(u71 + mm*(u72 + mm*(u73 + mm*(u73 + mm*(u72 + (u71 + mm*u70))))));
  uu8 = u80 + mm*(u81 + mm*(u82 + mm*(u83 + mm*(u84 + mm*(u84 + mm*(u83 + mm*(u82 + (u81 + mm*u80))))))));
  uu9 = u90 + mm*(u91 + mm*(u92 + mm*(u93 + mm*(u94 + mm*(u94 + mm*(u93 + mm*(u92 + (u91 + mm*u90))))))));
  
  series_sum = uu0 + yy*(uu1 + yy*(uu2 + yy*(uu3 + yy*(uu4 + yy*(uu5 + yy*(uu6 + yy*(uu7 + yy*(uu8 + ...
    yy*uu9))))))));
   
