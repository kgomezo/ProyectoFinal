set term png; set out 'Simulacion1.png'
set yr [-0.6: 0.6]
set xr [0: 60]
set grid
set xlabel 'Tiempo [s]' ; set ylabel '{/Symbol q} [rad]'
set title 'Angulo vs t sin fuerza aplicada'
b=sqrt(399)/20
A(x)=0.5*exp(-x/20)*cos(x*b)
o=sqrt(399)/798
B(x)=o*exp(-x/20)*sin(x*b)
C(x)=A(x)+B(x)
plot 'Simulacion1.txt' u 1:2 w l lt 4 t 'Runge-Kutta', C(x) t 'Teorico'

set term png; set out 'Simulacion21.png'
set grid
set xr [0: 65]
set xlabel 'Tiempo [s]'; set ylabel '{/Symbol q} [rad]'
set title 'Angulo vs t con F = 0.1'
plot 'Simulacion21.txt' u 1:2 w l lt 4 t 'F = 0.1'

set term png; set out 'Simulacion22.png'
set grid
set yr [-1: 1]
set xr [0: 65]
set xlabel 'Tiempo [s]'; set ylabel '{/Symbol q} [rad]'
set title 'Angulo vs t con F = 0.5'
plot 'Simulacion22.txt' u 1:2 w l lt 4 t 'F = 0.5'

set term png; set out 'Simulacion23.png'
set grid
set yr [-5: 5]
set xr [0: 65]
set xlabel 'Tiempo [s]'; set ylabel '{/Symbol q} [rad]'
set title 'Angulo vs t con F = 0.99'
plot 'Simulacion23.txt' u 1:2 w l lt 4 t 'F = 0.99'


set term png; set out 'Simulacion3.png'
set grid
set yr [-2: 2]
set xr [-5: 5]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title 'Sección de Poincare de Posicion vs Velocidad'
plot 'Simulacion3.txt' u 2:3 every 238::0 w p pt 7 t 'F = 1.2'

set term png; set out 'Simulacion611.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.4 [rad]'
plot'Simulacion611.txt' u 2:3 w p pt 7 t 'DT= 0.04'

set term png; set out 'Simulacion612.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.4 [rad]'
plot 'Simulacion612.txt' u 2:3 w p pt 7 t 'DT= 0.043'

set term png; set out 'Simulacion613.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.4 [rad]'
plot 'Simulacion613.txt' u 2:3 w p pt 7 t 'DT= 0.045'

set term png; set out 'Simulacion661.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.4 [rad]'
plot 'Simulacion611.txt' u 2:3 w p pt 7 t 'DT= 0.04', 'Simulacion612.txt' u 2:3 w p pt 7 t 'DT= 0.043', 'Simulacion613.txt' u 2:3 w p pt 7 t 'DT= 0.045'

set term png; set out 'Simulacion621.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.6 [rad]'
plot 'Simulacion621.txt' u 2:3 w p pt 7 t 'DT= 0.04'

set term png; set out 'Simulacion622.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.6 [rad]'
plot 'Simulacion622.txt' u 2:3 w p pt 7 t 'DT= 0.043'

set term png; set out 'Simulacion623.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.6 [rad]'
plot 'Simulacion623.txt' u 2:3 w p pt 7 t 'DT= 0.045'

set term png; set out 'Simulacion662.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.6 [rad]'
plot 'Simulacion621.txt' u 2:3 w p pt 7 t 'DT= 0.04', 'Simulacion622.txt' u 2:3 w p pt 7 t 'DT= 0.043', 'Simulacion623.txt' u 2:3 w p pt 7 t 'DT= 0.045'

set term png; set out 'Simulacion631.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.8 [rad]'
plot 'Simulacion631.txt' u 2:3 w p pt 7 t 'DT= 0.04'

set term png; set out 'Simulacion632.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.8 [rad]'
plot 'Simulacion632.txt' u 2:3 w p pt 7 t 'DT= 0.043'

set term png; set out 'Simulacion633.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.8 [rad]'
plot 'Simulacion633.txt' u 2:3 w p pt 7 t 'DT= 0.045'

set term png; set out 'Simulacion663.png'
set grid
set xr [-4: 4]
set yr [-3:3]
set xlabel '{/Symbol q} [rad]'; set ylabel '{/Symbol w} [rad/s]'
set title '{/Symbol D}{/Symbol q} vs t con {/Symbol q}_{0} = 0.8 [rad]'
plot 'Simulacion631.txt' u 2:3 w p pt 7 t 'DT= 0.04', 'Simulacion632.txt' u 2:3 w p pt 7 t 'DT= 0.043', 'Simulacion633.txt' u 2:3 w p pt 7 t 'DT= 0.045'
