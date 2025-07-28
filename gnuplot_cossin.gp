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

set xlabel "Angle $\\theta$"
set ylabel "Parameters"

set xrange [0:0.5*pi]
set yrange [-0.5:1]

set nokey

set output "Figures/gnuplot_cossq.tex"

f(x) = cos(x) * cos(x)
f2(x) = sin(x)**2
g(a,x) = cos(pi*a/10) * cos(x) * sin(x)

set title "Parameters $\\cos^2(\\theta)$, $\\sin^2(\\theta)$ and $\\sin\\cos(\\theta)\\sin(\\phi-\\phi_2)$"
plot f(x) with lines ls 1, f2(x) with lines ls 1, for [a=0:10] g(a,x) with lines ls 1 lt rgb "#FF00FF"


