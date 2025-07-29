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

set xlabel "Number of steps $n$"
set ylabel ""

set xrange [2:*]
set yrange [0.1:150]

set logscale y
set logscale x

set nokey

set term cairolatex pdf size 9cm,8cm color colortext font ",11"
set output "Figures/Classical_Variance.tex"

set title ""
plot for [corr=2:21] "Tabs/Correlated_Variance__ext.dat" using 1:corr with lines lw 2 lc rgb "#bbbbbbb" ,\
for [corr=2:21] "Tabs/Correlated_Variance_.dat" using 1:corr with lines lw 4 lc rgb monochrome_1(lin_int(corr,2,21,30,230))

f(c,x) = sqrt(1-sin(pi/2*c/19))*x
set ylabel "Variance $\\sigma^2$"
set term cairolatex pdf size 9.3cm,8cm color colortext font ",11"
set output "Figures/Quantum_Variance.tex"

set title ""
plot for [corr=0:19] f(corr,x) with lines lw 2 lc rgb "#bbbbbb" ,\
for [corr=2:21] "Tabs/Variances_Quantum_.dat" using 1:corr with lines lw 4 lc rgb monochrome_2(lin_int(corr,2,21,30,230))
