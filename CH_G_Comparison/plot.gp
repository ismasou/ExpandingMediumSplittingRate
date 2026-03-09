function $myplot(T,k) << EOD

    P = 16.0;
    CF = 4.0 / 3.0; CA = 3.0
    Pqg(z) = CF *(1.0 + (1.0 - z)**2) / z
    g = sqrt(4.0 * pi * 0.3);
    mDSqr = g**2 * T**2 * 1.5;
    # Harmonic Oscillator
    EulerGamma = 0.577216
    a = 4.0 * exp(-2.0 * EulerGamma + 1.0)
    muSqr = mDSqr
    qhat0 = g*2 * T * mDSqr / (2.0 * pi) * log(P / muSqr**0.5) / 1.8
    qhat1(QSqr, t, z) = 0.5 * ( 1.0 + (2.0 * CF / CA - 1.0)*z**2 + (1.0 - z)**2) * qhat0 * log(a*QSqr / muSqr)
    Qs(qhat, z) = sqrt(z * (1.0 - z) * P * qhat)

    # qhat2(t, z) = CA * qhat1(Qs(qhat1(1.0, t, z), z), t, z)
    qhat3(t, z) = qhat0 * ( CA + (2.0 * CF - CA) * z**2 + CA * (1.0 - z)**2) / 2.0
    I = sqrt(-1.0)
    w0(t, z) = sqrt(-I * qhat3(t, z) / (2.0 * P * z * (1.0 - z)))
    HO(t, z) = g**2 * Pqg(z) / (4.0 * pi**2 * P) * real( -w0(t,z) * tan(w0(t,z) * t))
    # qhat3(t, z) = qhat1(Qs(qhat2(t, z), z), t, z)
    # qhat4(t, z) = qhat1(Qs(qhat3(t, z), z), t, z)
    # qhat4(t, z) = qhat3(Qs(qhat3(t, z), z), z)
    # print qhat1(1.0, 0.2, 0.5), qhat2(0.2, 0.5), qhat3(0.2, 0.5), qhat4(0.2, 0.5)

    unset log x
    AMY= 0.0194
    z = k / 16.0;
    TT=sprintf("%.1f_%.0f",T,k)
    set title sprintf("T=%.1f, k=%.0f",T,k)
    factor= Pqg(z) / P
    if (T == 0.2 && k == 3.0) {
        set xrange[0:5]
        set yrange[0:0.02]
        AMY = 0.0077
    }
    if (T == 0.4 && k == 3.0) {
        set xrange[0:3]
        set yrange[0:0.04]
        AMY = 0.0194
    }
    if (T == 0.2 && k == 8.0) {
        set xrange[0:5]
        set yrange[0:0.003]
        AMY = 0.0194
    }
    if (T == 0.4 && k == 8.0) {
        set xrange[0:5]
        set yrange[0:0.01]
        AMY = 0.0194
    }
    set key top left
    dt(x) = 2.0 * P * z * (1.0 - z) / mDSqr * x;
    lambda = (g**2*T**2) * (2.0 * P *z * (1.0 - z)) / mDSqr;
    InvGevTofm = 0.197326 ;
    print T,k
    plot\
        "./CH_Gale-".TT.".txt" u 1:2 w l ls 1 lw 1 ti "Caron-Huot Gale",\
        "./rate-".TT.".dat" u (InvGevTofm * dt($1)):(g**4 * T * Pqg(z) / pi * $2 / P) w l ls 5 lw 1 dt 2 ti "SplittingRates",\
        "./rate-2d-".TT.".dat" u (InvGevTofm * dt($1)):(g**4 * T * Pqg(z) / pi * $2 / P) w l ls 6 lw 1 dt 4 ti "2D Integral",\
        "./rate-opacity-".TT.".dat" u (InvGevTofm * dt($1)):(g**4 * T * Pqg(z) / pi * $2 / P) w l ls 7 lw 1 ti "opacity",\
        HO(x / InvGevTofm, z) w l ls 3 lw 0.6 dt (2,2) ti "HO"
    return 1
EOD

set term pngcairo size 1920,1080 enhanced font ",20" lw 8
set output "./CG_G_Comparison.png"
set multiplot layout 2,2
T = 0.2; k = 3;
x = $myplot(T,k)
T = 0.4; k = 3;
x = $myplot(T,k)
T = 0.2; k = 8;
x = $myplot(T,k)
T = 0.4; k = 8;
x = $myplot(T,k)
unset multiplot
set output
