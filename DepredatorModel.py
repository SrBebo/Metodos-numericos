import numpy as np

def rungekutta2_fg(f,g,t0,x0,y0,h,muestras):
    tamano = muestras +1
    tabla = np.zeros(shape=(tamano,3),dtype=float)
    tabla[0] = [t0,x0,y0]
    ti = t0
    xi = x0
    yi = y0
    for i in range(1,tamano,1):
        K1x = h * f(ti,xi,yi)
        K1y = h * g(ti,xi,yi)
        
        K2x = h * f(ti+h, xi + K1x, yi+K1y)
        K2y = h * g(ti+h, xi + K1x, yi+K1y)

        xi = xi + (1/2)*(K1x+K2x)
        yi = yi + (1/2)*(K1y+K2y)
        ti = ti + 2*h
        tabla[i] = [ti,xi,yi]
    tabla = np.array(tabla)
    return(tabla)

# PROGRAMA ------------------

# INGRESO
# Parámetros de las ecuaciones
a = 0.3#0.5
b = 0.2#0.7
c = 0.5#0.35
d = 0.5#0.35

# Ecuaciones
f = lambda t,x,y : a*x -b*x*y
g = lambda t,x,y : -c*y + d*x*y

# Condiciones iniciales
t0 = 0
x0 = 2
y0 = 1

# parámetros del algoritmo
h = 0.5
muestras = 102


# PROCEDIMIENTO
tabla = rungekutta2_fg(f,g,t0,x0,y0,h,muestras)
ti = tabla[:,0]
xi = tabla[:,1]
yi = tabla[:,2]

# SALIDA
np.set_printoptions(precision=6)
print(' [ ti, xi, yi]')
print(tabla)

# Grafica tiempos vs población
import matplotlib.pyplot as plt

plt.plot(ti,xi, label='xi presa')
plt.plot(ti,yi, label='yi predador')

plt.title('Modelo predador-presa')
plt.xlabel('t tiempo')
plt.ylabel('población')
plt.legend()
plt.grid()
plt.show()

# gráfica xi vs yi
plt.plot(xi,yi)

plt.title('Modelo presa-predador [xi,yi]')
plt.xlabel('x presa')
plt.ylabel('y predador')
plt.grid()
plt.show()
