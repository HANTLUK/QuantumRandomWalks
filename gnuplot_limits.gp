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
set xtics textcolor rgb "#808080"
set ytics textcolor rgb "#808080"
set y2tics textcolor rgb "#808080"

rgb(r,g,b) = int(r*65536 + g*256 + b)
lin_int(value,min,max,val_min,val_max) = val_min + (val_max-val_min)/(max-min)*(value-min)
monochrome_1(b) = rgb(b,100,100)
monochrome_2(b) = rgb(100,b,100)
monochrome_3(b) = rgb(100,100,b)

set xlabel "$x$"
set ylabel "Distribution $f(x)$"

set yrange [0:*]

set nokey

set output "Figures/Classical_Distribution.tex"

set title "Limiting Distribution of the classically correlated random walk"
plot for [corr=2:10] "Tabs/Classical_Distribution.dat" using 1:corr with lines lw 4 lc rgb monochrome_1(lin_int(corr,2,10,100,255))


set xrange [-1:1]
set yrange [0:*]

set nokey

ang(par) = pi/2*par/19
f(par,x) = (sqrt(1-cos(ang(par))**2))/(pi*(1-x**2)*sqrt(cos(ang(par))**2-x**2))

set output "Figures/Quantum_Distribution_Symmetric.tex"

set title "Limiting Distribution of the symmetric quantum random walk"
plot for [corr=1:18] [-cos(ang(corr)):cos(ang(corr))] f(corr,x) with lines lw 4 lc rgb monochrome_2(lin_int(corr,1,18,100,255))


set xrange [-0.4:0.4]
set yrange [10**(-1):10**2]
set logscale y

set nokey

set output "Figures/Quantum_Distribution.tex"

set title "Limiting Distribution of the quantum random walk"
plot for [corr=2:21] "Tabs/Quantum_Distribution.dat" using 1:corr with lines lw 4 lc rgb monochrome_3(lin_int(corr,2,21,100,255))
