# GNUPLOT TEMPLATE FILE(c) - Federico Moretta
#!/usr/bin/env/gnuplot

# TAKE FILE NAME AND PATH FROM COMMAND LINE
# to take a file from the command line, use the following command:
# $ gnuplot -e "filename='filename.dat'" PLOT.plt
#reset
# if ( !exists("filename") ) quit("No file name given")

# SET TERMINAL TYPE FOR VECTORIAL IMAGES
# set terminal postscript enhanced size 1000,700
# set output "reg_res".'.eps'

set terminal svg enhanced size 1400,600 font "Verdana,12"
set encoding 'utf8'
set output "reg_res_2".'.svg'

# to save it in png, switch these two terminals
# set terminal pngcairo enhanced size 1000,1000
# set output 'gnuplot_figure.png'

# PLOTTING PARAMETERS
set autoscale                                   # scale axes automatically
set size square                                 # make the plot square

unset log                                       # remove any log-scaling
unset label                                     # remove any previous labels

set xtic auto                                   # set xtics automatically
set ytic auto                                   # set ytics automatically
set grid                                        # ENABLE GRID

colour = 2
point = 2

set multiplot layout 2,5

set key left top                    
set xlabel 'X_T'
set ylabel 'D(X_T°-X_T)'
set title 'k_h regression'
set xtic auto
plot "kinetics_hyd.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_hyd.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel 'D(S_1° - S_1)+k_hX_T'
set ylabel '{/Symbol m}_1X_1' enhanced
set title 'v_1 regression'
set xtic 0,0.06,0.35
plot "kinetics_k1.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k1.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel '{/Symbol m}_2' enhanced
set ylabel 'qCH_4/X_2'
set title 'v_6 regression'
set xtic 0,0.04,0.23
plot "kinetics_k6.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k6.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel 'SCO_2-K_HP_C'
set ylabel 'qCO_2'
set title 'k_La regression'
set xtic -2.4, 0.3, -1
plot "kinetics_kLa.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_kLa.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set ylabel 'qCO_2-D(C°-C)'
set xlabel 'D(S_1°-S_1)+k_hX_T'
set title 'v4/v1 regression'
set xtic auto
plot "kinetics_k4k1_k5k6_1.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k4k1_k5k6_1.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set ylabel 'qCO_2-D(C°-C)'
set xlabel 'qM'
set title 'v5/v6 regression'
set xtic auto
plot "kinetics_k4k1_k5k6_2.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k4k1_k5k6_2.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left bottom
set xlabel 'D(S_2°-S_2)'
set ylabel 'qM'
set title 'v6/v3 regression'
set xtic -18, 3, 0
plot "kinetics_k6k3_k2k1_1.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k6k3_k2k1_1.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel 'D(S_1°-S_1)+ k_hX_T'
set ylabel 'qM'
set title 'v6v2/v3v1 regression'
set xtic auto
plot "kinetics_k6k3_k2k1_2.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_k6k3_k2k1_2.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel 'DS_1'
set ylabel 'S_1'
set title '{/Symbol m}m_1, Ks_1, Cd_1 regression' enhanced
set xtic 0, 0.035, 0.2
plot "kinetics_S1_1.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_S1_1.txt" using 1:3 title "mod" \
w lp pt point lc colour

set key left top
set xlabel 'D'
set ylabel 'S_1'
set title '{/Symbol m}m_1, Ks_1, Cd_1 regression' enhanced
set xtic 0, 0.035, 0.2
plot "kinetics_S1_2.txt" using 1:2 title "exp" \
w point pt 6 lc rgb "black", "kinetics_S1_2.txt" using 1:3 title "mod" \
w lp pt point lc colour

unset multiplot