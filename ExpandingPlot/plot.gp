set terminal cairolatex eps color font ",20in" size 6,11 standalone lw 4 ps 0.15 background "white" header '\
\newcommand{\hl}[1]{\setlength{\fboxsep}{0.75pt}\colorbox{white}{#1}}\
\usepackage{amsmath}\
\newcommand{\mysize}{\small}\
'
set size square

load "./common.gp"

OutputFold= "Temp"
system("mkdir -p ".OutputFold)

fname = "InitialTime"
File = OutputFold."/".fname
set output File.".tex"

set lmargin 9
set ylabel '\mysize Radiation rate: $\frac{d\Gamma}{dz}(P,z)$' offset 2,0
set xlabel '\mysize Time $[$fm$/$c$]$' offset 0,0.5
set key top left offset 0,-0.1 reverse Left maxrows 5 width -5 spacing 0.8
Na = 3; Nt = 2;
t0 = 1e-4;
set format x '\mysize $%g$'
set format y '\mysize $%g$'
set yrange [0:2.1]
set xrange [0:0.1]
set multiplot layout 2,1

do for [zval in "0.01"]{
    print zval
    set yrange[0:0.076*Pgg(zval)]
    set xrange[0:(InvGevTofm * dt(15.55, zval))]
    eveFactor = 2
    set label 3 at graph 0.70 , 0.72 '\mysize $z = '.zval.'$' front
    plot\
        keyentry w l ls 10 lw 2 ti '\mysize Brick~~~~~~~',\
        for[i=Nt:1:-1] keyentry w p ls 1 lw 2 pt PoinType[i] ps 10 ti gprintf('\mysize $t_0 =%g$ fm/c', t0Vals[i]),\
        for[i=Nt:1:-1] keyentry w l ls 1 dt Nt-i+1 ti gprintf('\mysize OE $t_0 =%g$ fm/c', t0Vals[i]),\
        "./Full/rate_P".P."_z".zval."_T".T.".dat"\
        u @USEtime:@USEX w l ls 10 lw 2 ti '',\
        for[k=3:Na] keyentry w l ls 1+k ti gprintf('\mysize $\alpha =%g$', alpha[k]),\
        for[k=3:Na] for[i=1:Nt] sprintf("./Opacity/bjorken_P".P."_z".zval."_T".T."_alpha%g_t0%g.dat", alpha[k], t0Vals[i])\
            u @USE eve eveFactor w l ls 1+k dt Nt-i+1 ti "",\
        for[k=3:Na] for[i=1:Nt] sprintf("./Expanding/rate_P".P."_z".zval."_T".T."_alpha%g_t0%g.dat", alpha[k], t0Vals[i])\
            u (@USEtime < 0.1 ? @USEtime:1/0):@USEX w p ls 1+k lw 2 pt PoinType[i] ps 10 ti "",\
        for[k=3:Na] for[i=1:Nt] sprintf("./Expanding/rate_P".P."_z".zval."_T".T."_alpha%g_t0%g.dat", alpha[k], t0Vals[i])\
            u @USE eve eveFactor w p ls 1+k lw 2 pt PoinType[i] ps 10 ti "",\
            keyentry ti " ",\
            keyentry ti " ",\
            keyentry ti " "
#
}

do for [zval in "0.5"]{
    print zval
    set yrange[0:0.03*Pgg(zval)]
    set xrange[0:(InvGevTofm * dt(5.55, zval))]
    eveFactor = 1
    set label 3 at graph 0.70 , 0.62 '\mysize $z = '.zval.'$' front
    plot\
        for[i=Nt:1:-1] keyentry w p ls 1 lw 2 pt 2+(i-1)*4 ps 10 ti gprintf('\mysize $t_0 =%g$ fm/c', t0Vals[i]),\
        for[i=Nt:1:-1] keyentry w l ls 1 dt Nt-i+1 ti gprintf('\mysize HO $t_0 =%g$ fm/c', t0Vals[i]),\
        "./Full/rate_P".P."_z".zval."_T".T.".dat"\
        u @USE w l ls 10 lw 2 ti '',\
        for[k=3:Na] keyentry w l ls 1+k ti gprintf('\mysize $\alpha = %g$', alpha[k]),\
        for[k=3:Na] for[i=1:Nt] sprintf("./HO/rate_P".P."_z".zval."_T".T."_alpha%g_t0%g.dat", alpha[k], t0Vals[i])\
            u ($1+t0Vals[i]):($2) w l ls 1+k dt Nt-i+1 ti "",\
        for[k=3:Na] for[i=1:Nt] sprintf("./Expanding/rate_P".P."_z".zval."_T".T."_alpha%g_t0%g.dat", alpha[k], t0Vals[i])\
            u @USE w p ls 1+k lw 2 pt PoinType[i] ps 10 ti "",\
        keyentry w l ls 10 lw 2 ti '\mysize Brick~~~~~~~',\
        keyentry ti " ",\
        keyentry ti " "
#
}

unset multiplot
set output

sys("latexmk -lualatex -f -ps -jobname=".File." ".File.".tex " )
sys("ps2pdf ".File.".ps ".File.".pdf")
sys("mv ".File.".pdf ".fname.".pdf")
sys("rm -rf Temp")
# sys("clear")

