set term png; set out 'Simulacion4.png'
set grid
set logscale y
set xlabel 'Tiempo [s]'; set ylabel '{/Symbol D}{/Symbol q} [rad]'
set title '{/Symbol D}{/Symbol q} vs t con F = 1.2'
A(x)=c*x+b
fit A(x) 'Simulacion4.txt' u 1:2 via c, b
A(x) = exp(c*x+b)
fit A(x) 'Simulacion4.txt' u 1:2 via c, b
set yrange [0.00001:100]
set xr [0:65]
plot 'Simulacion4.txt' u 1:2 w l lt 4 t 'F = 1.2', A(x) t 'e^{{/Symbol l}x}'

set term png; set out 'Simulacion5.png'
set grid
set xr [0:65]
set logscale y
set xlabel 'Tiempo [s]'; set ylabel '{/Symbol D}{/Symbol q} [rad]'
set title '{/Symbol D}{/Symbol q} vs t con q variable'
A(x)=c*x+b
fit A(x) 'Simulacion5.txt' u 1:2 via c, b
A(x) = exp(c*x+b)
fit A(x) 'Simulacion5.txt' u 1:2 via c, b
plot 'Simulacion5.txt' u 1:2 w l lt 4 t 'F = 1.2', A(x) t 'e^{{/Symbol l}x}'
