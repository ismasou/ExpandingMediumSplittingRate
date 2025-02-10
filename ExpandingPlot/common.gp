P="16"; T="0.4";

CF = 4.0 / 3.0; CA = 3.0
Pgg(z) = 2.0 * CA * (1.0 - z * (1.0 - z))**2 / (z * (1.0 - z))
g = sqrt(4.0 * pi * 0.3);
mDSqr = g**2 * T**2 * 1.5;

#######################
# HARMONIC OSCILLATOR #
#######################
# EulerGamma = 0.577216
# a = 4.0 * exp(-2.0 * EulerGamma + 1.0)
# muSqr = mDSqr
# qhat0 = g**2 * T * mDSqr / (2.0 * pi) * log(P / muSqr**0.5) / 1.8
# qhat1(QSqr, t, z) = 0.5 * ( 1.0 + (2.0 * CF / CA - 1.0)*z**2 + (1.0 - z)**2) * qhat0 * log(a*QSqr / muSqr)
# Qs(qhat, z) = sqrt(z * (1.0 - z) * P * qhat)
#
# # qhat2(t, z) = CA * qhat1(Qs(qhat1(1.0, t, z), z), t, z)
# qhat3(t, z) = qhat0 * ( CA + CA * z**2 + CA * (1.0 - z)**2) / 2.0
# I = sqrt(-1.0)
# w0(t, z) = sqrt(-I * qhat3(t, z) / (2.0 * P * z * (1.0 - z)))
# HO(t, z) = g**2 * Pgg(z) / (4.0 * pi**2 * P) * real( -w0(t,z) * tan(w0(t,z) * t))
# qhat3(t, z) = qhat1(Qs(qhat2(t, z), z), t, z)
# qhat4(t, z) = qhat1(Qs(qhat3(t, z), z), t, z)
# qhat4(t, z) = qhat3(Qs(qhat3(t, z), z), z)
# print qhat1(1.0, 0.2, 0.5), qhat2(0.2, 0.5), qhat3(0.2, 0.5), qhat4(0.2, 0.5)


#############
# FUNCTIONS #
#############
dt(x,z) = 2.0 * P * z * (1.0 - z) / mDSqr * x;
lambda(z) = (g**2*T**2) * (2.0 * P *z * (1.0 - z)) / mDSqr;
InvGevTofm = 0.197326 ;

array alpha=[1e-2, 4e-1, 1e0, 1e0] # 1e-3
# array t0Vals = [0.01, 0.1, 0.3]
array t0Vals = [0.1, 0.3, 0.3]
Temp(t,a,t0) = ((t0 / (t + t0 + 10.0))**(a/ 3.0))
array PoinType = [6, 2]


USEtime="(InvGevTofm * dt($1, zval))"
USEX="(g**4 * T * Pgg(zval) / pi * $2)"
# USEX2="(g**4 * T * Temp(dt($1, zval), alpha[k], @t0) * Pgg(zval) / pi * $2)"
USE=USEtime.":".USEX
# USE2=USEtime.":".USEX2
USEX2="(g**4 * T * Pgg(zval) / pi * $2 / Temp(InvGevTofm * dt($1, zval), alpha[k], t0Vals[i]))"
USEX2="(g**4 * T * Pgg(zval) / pi * $2 * Temp(InvGevTofm * dt($1, zval), alpha[3], t0Vals[2]))"
USE2=USEtime.":".USEX2

posi = 4
set style line 1 ps posi lc rgb '#110141';   set style line 2 ps posi lc rgb '#710162';   set style line 11 ps posi lc rgb '#a12a5e';
set style line 5 ps posi lc rgb '#ef6a32';   set style line 12 ps posi lc rgb '#fbbf45';
set style line 7 ps posi lc rgb '#aad962';   set style line 9 ps posi lc rgb '#017351';
set style line 3 ps posi lc rgb '#01545a';   set style line 10 ps posi lc rgb '#26294a';   set style line 8 ps posi lc rgb '#1a1334';

set style line 9 ps posi lt 1 lc rgb '#7FC97F' # pale green
set style line 9 ps posi lt 1 lc rgb '#BEAED4' # pale purple
set style line 9 ps posi lt 1 lc rgb '#FDC086' # pale orange
set style line 9 ps posi lt 1 lc rgb '#FFFF99' # pale yellow
set style line 9 ps posi lt 1 lc rgb '#F0027F' # magenta
set style line 6 ps posi lt 1 lc rgb '#BF5B17' # brown
set style line 9 ps posi lt 1 lc rgb '#666666' # grey

set style line 2 lw 2 ps posi lc rgb '#386CB0' # blue
set style line 3 lw 2 ps posi lc rgb '#03c383';
set style line 4 lw 2 ps posi lc rgb '#ed0345';
set style line 5 lw 2 ps posi lc rgb '#fbbf45';
