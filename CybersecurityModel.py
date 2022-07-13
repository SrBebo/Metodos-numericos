
import numpy as np

def rungekutta2_fgw(f,g,w,t0,x0,y0,z0,h,muestras):
    tamano = muestras +1
    tabla = np.zeros(shape=(tamano,4),dtype=float)
    tabla[0] = [t0,x0,y0,z0]
    ti = t0
    xi = x0
    yi = y0
    zi = z0

    for i in range(1,tamano,1):
        K1x = h * f(ti,xi,yi,zi)
        K1y = h * g(ti,xi,yi,zi)
        K1z = h * w(ti,xi,yi,zi)
        
        K2x = h * f(ti+h/2, xi + K1x/2, yi+K1y/2, zi +K1z/2)
        K2y = h * g(ti+h/2, xi + K1x/2, yi+K1y/2, zi +K1z/2)
        K2z = h * w(ti+h/2, xi + K1x/2, yi+K1y/2, zi +K1z/2)

        K3x = h * f(ti+h/2, xi + K2x/2, yi+K2y/2, zi +K2z/2)
        K3y = h * g(ti+h/2, xi + K2x/2, yi+K2y/2, zi +K2z/2)
        K3z = h * w(ti+h/2, xi + K2x/2, yi+K2y/2, zi +K2z/2)

        K4x = h * f(ti+h, xi + K3x, yi+K3y, zi +K3z)
        K4y = h * g(ti+h, xi + K3x, yi+K3y, zi +K3z)
        K4z = h * w(ti+h, xi + K3x, yi+K3y, zi +K3z)
        
        xi = xi + (1/6)*(K1x+2*K2x+2*K3x+K4x)
        yi = yi + (1/6)*(K1y+2*K2y+2*K3y+K4y)
        zi = zi + (1/6)*(K1z+2*K2z+2*K3z+K4z)
        ti = ti + h
        
        tabla[i] = [ti,xi,yi,zi]
    tabla = np.array(tabla)
    return(tabla)
  
# Programa
# Parámetros de las ecuaciones

binfecta = 0.5
grecupera = 0.01#[0.00001-0.009]
# Ecuaciones
f = lambda t,S,I,R : -binfecta*S*I#Susceptibilidad de equipos
g = lambda t,S,I,R : binfecta*S*I - grecupera*I#Cantidad de dispositivos infectados
w = lambda t,S,I,R : grecupera*I#Dispositivos recuperados
# Condiciones iniciales
t0=1
S0 = 1.0#tasa de susceptibilidad a un posible ataque
I0 = 0.01#tasa de infectados [0.001-0.1]
R0 = 0.1#tasa de recuperacion de infeccion[0.00001-0.01] 
#Paramentros del algoritmo
h = 1#numero de dias 
muestras = 100
# PROCEDIMIENTO
tabla = rungekutta2_fgw(f,g,w,t0,S0,I0,R0,h,muestras)
ti = tabla[:,0]
Si = tabla[:,1]
Ii = tabla[:,2]
Ri = tabla[:,3]
# SALIDA
np.set_printoptions(precision=6)
print(' [ ti, Si, Ii, Ri]')
print(tabla)

# Grafica tiempos vs equipos infectados
import matplotlib.pyplot as plt
plt.plot(ti,Si, label='Susceptibles')
plt.plot(ti,Ii, label='Infectados')
plt.plot(ti,Ri, label='Recuperados')
plt.title('Modelo SIR')
plt.xlabel('t tiempo')
plt.ylabel('equipos')
plt.legend()
plt.grid()
plt.show()
