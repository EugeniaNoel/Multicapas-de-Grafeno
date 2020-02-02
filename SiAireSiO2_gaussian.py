import numpy as np
from matplotlib import pyplot as plt

#Cargo las conductividades que calcule con los programas 'conductividad.py' y 'conductividad_0.py'
v = np.loadtxt('sigma_xx 5T adim large.txt')[0]
 
#constantes
c = 2.99792458*10**10 #cm/s
e0 = 8.8541878176*10**(-12) #N/V2
u0 = 4*np.pi*10**(-7) #Tm/A
h = 6.626070*10**(-34) #kg m2/s
h_barra = h/(2*np.pi)
e = 1.6*10**(-19) #C

#parametros 
ep1= 12 #Si
n1 = np.sqrt(ep1)
ep2= 1 #Aire
n2 = np.sqrt(ep2)
ep3 = 4 #SiO2
n3 = np.sqrt(ep3)

tita1 = 60*np.pi/180 # angulo de incidencia en radianes
cos_tita1 = np.cos(tita1)
sen_tita1 = np.sin(tita1) 
cos_tita2 = np.sqrt(1-(ep1/ep2)*sen_tita1**2, dtype=complex) # el coseno del angulo transmitido en 2 va a ser complejo
cos_tita3 = np.sqrt(1-(ep1/ep3)*sen_tita1**2, dtype=complex)

def k(omega, ep):   # para calcular el k en cada medio
    k = np.sqrt(ep)*omega/c
    return(k)

omega = 2*np.pi*np.array(v)

k1 = k(omega, ep1)
q = k1*sen_tita1
k1z = k1*cos_tita1

#k0 = omega/c

k2 = k(omega, ep2) 
k2z = k2* cos_tita2 
#k2z = 1j*k0*np.sqrt((n1**2)*sen_tita1**2-n2**2) #lo sque de unas cuentas de wikipedia https://en.wikipedia.org/wiki/Total_internal_reflection#Evanescent_wave_(qualitative_explanation)
#Da lo mismo que el k2z que deje
#k2z1 = k2*cos_tita2
k3 = k(omega, ep3)
k3z = k3*cos_tita3

# =============================================================================
# #Grafico para ver si Re(q)>k y que im(kz) > 0
# plt.figure()
# plt.plot(omega, k2, 'o', label='k2')
# plt.plot(omega, q, 'o', label='q')
# plt.plot(omega, np.imag(k2z), 'o', label='im(k2z)')
# plt.legend(loc='best')
# plt.xlabel('omega')
# plt.ylabel('numero de onda')
# 
# =============================================================================
#alfa2 = ((ep2*ep3)/(ep2+ep3))*(omega/c)**2 # chequee que sea igual a q3**2 pq esto sale de la condicion de contorno

#Matriz de desfasaje en un medio
def desfasaje(Dz, kz): #Esta en funcion de la distancia que se recorre en ese medio y el numero de onda angular de ese medio
    D = np.zeros((4,4), dtype= np.complex)
    ida = np.exp(-1j*kz*Dz)
    vuelta = np.exp(1j*kz*Dz)
    D[0,0] = ida    
    D[2,2] = ida
    D[1,1] = vuelta
    D[3,3] = vuelta
    return(D)
    

#%% Vector con la matriz de transición de medio 1 al 2 en funcion de la frecuencia
# Esribo la matriz de transfrencia considerando conductividades igual 0

cte5= (n2/n1)*(cos_tita2/cos_tita1)   
cte1= cte5*ep1/ep2
T = []
for i in range(len(omega)):
    T_i= np.array([[ 1+cte1 , 1-cte1, 0 ,  0 ], 
       [ 1-cte1 , 1+cte1 , 0 , 0 ],
       [ 0 , 0 , 1+cte5 , 1-cte5],
       [ 0 , 0 , 1-cte5 , 1+cte5 ]], dtype = complex)
    T_i=np.multiply(0.5,T_i, dtype=complex)
    T.append(T_i)
  
#%%  Vector con la matriz de transición de medio 2 al 3 en funcion de la frecuencia
# Como tengo la capa de grafeno en el medio, uso la conductividad del grafeno 
    
Rxx = np.array(np.loadtxt('sigma_xx 5T adim large.txt')[1]) 
Ixx = np.array(np.loadtxt('sigma_xx 5T adim large.txt')[2])
#Rxy = np.array([0]*len(v-1))  # caso B=0
#Ixy = np.array([0]*len(v-1))  # caso B=0
Rxy = np.array(np.loadtxt('sigma_xy 5T adim large.txt')[1])
Ixy = np.array(np.loadtxt('sigma_xy 5T adim large.txt')[2])


#Ya reemplazando que sigma_xx = sigma_yy y sigma_xy = -sigma_yx
lv = len(v)
E = np.zeros(lv, dtype=object)
for i in range(lv):
    sigma_xx = Rxx[i]+Ixx[i]*1j 
    sigma_xy = Rxy[i]+Ixy[i]*1j
    cte5= k3z[i]/k2z[i]
    cte1= cte5*(ep2/ep3)
    cte3 = 4*np.pi/137.036
    cte2= cte3*k3z[i]*c/(ep3*omega[i])
    cte4= cte3*cte5/ep3
    cte6= omega[i]*cte3/(k2z[i]*c)
    E_i= np.array([[ (1+cte1+sigma_xx*cte2), (1-cte1-sigma_xx*cte2), (sigma_xy*cte3) ,  (sigma_xy*cte3) ], 
           [ (1-cte1+sigma_xx*cte2) , (1+cte1-sigma_xx*cte2), (sigma_xy*cte3), (sigma_xy*cte3) ],
           [ -(sigma_xy*cte4) , (sigma_xy*cte4) , (1+cte5+sigma_xx*cte6) , (1-cte5+sigma_xx*cte6) ],
           [ (sigma_xy*cte4) , -(sigma_xy*cte4) , (1-cte5-sigma_xx*cte6) , (1+cte5-sigma_xx*cte6) ]], dtype = complex)
    E_i = np.multiply(0.5, E_i, dtype=complex)
    E[i] = E_i


d = 5*10**-4     #Ancho de la capa de Aire en cm
Mfinal = []      #Matriz de trasnmicion final de la configuracion SiAireSiO2

for i in range(lv):
    D = desfasaje(d, k2z[i])
    D1 = np.dot(T[i], D)
    M0 = np.dot(D1, E[i])
    Mfinal.append(M0)

#%% Calculo las reflectancias y transmitancias en funcion de la frecuencia
    
aN_hy = []
cN_hy = []
b1_hy = []
d1_hy = []
aN_ey = []
cN_ey = []
b1_ey = []
d1_ey = []

for n in range(len(Mfinal)):
    M = Mfinal[n]
    for i in range(len(M[:,0])):
        for j in range(len(M[0,:])):
            globals()['M' +str(i+1)+ str(j+1)] = M[i,j]
    #Para incidente TM
    aN_hyn= -1*M33/(M13*M31-M11*M33)
    cN_hyn = M31 /(M13*M31-M11*M33)
    b1_hyn = (M31*M23-M21*M33)/(M13*M31-M11*M33)
    d1_hyn = (M43*M31-M41*M33)/(M13*M31-M11*M33)
    aN_hy.append(aN_hyn)
    cN_hy.append(cN_hyn)
    b1_hy.append(b1_hyn)
    d1_hy.append(d1_hyn)
    #Para incidente TE
    aN_eyn= -1*M13/(M33*M11-M13*M31)
    cN_eyn = M11 /(M33*M11-M31*M13)
    b1_eyn = (M11*M23-M21*M13)/(M33*M11-M13*M31)
    d1_eyn = (M43*M11-M13*M41)/(M33*M11-M31*M13) 
    aN_ey.append(aN_eyn)
    cN_ey.append(cN_eyn)
    b1_ey.append(b1_eyn)
    d1_ey.append(d1_eyn)

a = cos_tita3/cos_tita1


RefTM_hy = (np.abs(b1_hy)**2)
RefTE_hy = (np.abs(d1_hy)**2)*ep1
TransTM_hy = a*(np.abs(aN_hy)**2)*(n1/n3)
TransTE_hy = a*(np.abs(cN_hy)**2)*(n3*n1)

RefTM_ey = (np.abs(b1_ey)**2)/ep1
RefTE_ey = (np.abs(d1_ey)**2)
TransTM_ey = a*(np.abs(aN_ey)**2)/(n3*n1)
TransTE_ey = a*(np.abs(cN_ey)**2)*(n3/n1)

font = {'family': 'arial', 'weight': 'normal','size': 14 }

plt.figure()
plt.grid()
plt.plot(v,RefTM_hy, 'r-')
plt.plot(v,TransTM_hy, 'b-' )
plt.plot(v,RefTE_hy, 'g-')
plt.plot(v,TransTE_hy, 'y-' )
plt.legend(('R TM', 'T TM' , 'R TE', 'T TE' ),loc='best', fontsize=12 )
plt.title('Incidencia TM', fontdict=font)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.ylim((-0.1,1.1))
plt.xlabel(r'$\nu$ (Hz)', fontdict=font)
plt.ylabel('Fracción de la intensidad incidente', fontdict=font)

plt.figure()
plt.grid()
plt.plot(v,RefTM_ey, 'r-')
plt.plot(v,TransTM_ey, 'b-' )
plt.plot(v,RefTE_ey, 'g-')
plt.plot(v,TransTE_ey, 'y-' )
plt.legend(('R TM', 'T TM', 'R TE', 'T TE' ),loc='best', fontsize=12 )
font = {'family': 'arial', 'weight': 'normal','size': 14 }
plt.title('Incidencia TE', fontdict=font)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.ylim((-0.1,1.1))
plt.xlabel(r'$\nu$ (Hz)', fontdict=font)
plt.ylabel('Fracción de la intensidad incidente', fontdict=font)

np.savetxt('5T TM adim large.txt', (v, RefTM_hy, TransTM_hy, RefTE_hy, TransTE_hy))
np.savetxt('5T TE adim large.txt', (v, RefTM_ey, TransTM_ey, RefTE_ey, TransTE_ey))
