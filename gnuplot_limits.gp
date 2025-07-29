set term cairolatex pdf size 8cm,8cm color colortext font ",11"

set decimalsign "." # für den input
#               Style
# !!!___________________________ !!!
set grid xtics
set grid ytics
# set grid mxtics
# set grid mytics
set style line 80 linetype 1 linecolor rgb "#888888"
set style line 81 linetype 1 linecolor rgb "#808080" linewidth 0.5
set border back linestyle 80
set grid back linestyle 81
set xtics textcolor rgb "#000000"
set ytics textcolor rgb "#000000"
set y2tics textcolor rgb "#000000"

rgb(r,g,b) = int(r*65536 + g*256 + b)
lin_int(value,min,max,val_min,val_max) = val_min + (val_max-val_min)/(max-min)*(value-min)
monochrome_1(b) = rgb(b,230-b,230-b)
monochrome_2(b) = rgb(b,230-b,b)
monochrome_3(b) = rgb(230-b,230-b,b)

set xlabel "$x$"
set ylabel ""
set yrange [0:*]
set nokey
set term cairolatex pdf size 9cm,8cm color colortext font ",11"
set output "Figures/Classical_Distribution.tex"
set title ""
plot for [corr=2:10] "Tabs/Classical_Distribution.dat" using 1:corr with lines lw 4 lc rgb monochrome_1(lin_int(corr,2,10,30,230))


n = 7
ang(par) = pi/2*par/(n+1)
f(par,x) = (sqrt(1-cos(ang(par))**2))/(pi*(1-x**2)*sqrt(cos(ang(par))**2-x**2))
set ylabel "Distribution $f(x)$"
set term cairolatex pdf size 9.3cm,8cm color colortext font ",11"
set output "Figures/Quantum_Distribution_Symmetric_Cut.tex"
set xrange[-1:1]
part = 0.993
do for [corr=1:n] {
    set arrow from -part*cos(ang(corr)),0.01 to -part*cos(ang(corr)),f(corr,-part*cos(ang(corr))) nohead lw 4 lc rgb monochrome_2(lin_int(corr,1,n,30,230))
    set arrow from part*cos(ang(corr)),0.01 to part*cos(ang(corr)),f(corr,part*cos(ang(corr))) nohead lw 4 lc rgb monochrome_2(lin_int(corr,1,n,30,230))
}
set title ""
set samples 1000
plot for [corr=1:n] [-cos(ang(corr)):cos(ang(corr))] f(corr,x) with lines lw 4 lc rgb monochrome_2(lin_int(corr,1,n,30,230))
unset arrow

n = 21
set xrange [-0.31:0.31]
set yrange [0.05:50]
set logscale y
set nokey
set term cairolatex pdf size 9.3cm,8cm color colortext font ",11"
set output "Figures/Quantum_Distribution.tex"
set title ""
plot for [corr=2:n] "Tabs/Quantum_Distribution.dat" using 1:corr with lines lw 4 lc rgb monochrome_3(lin_int(corr,2,n,30,230))
