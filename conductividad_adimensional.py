import numpy as np
import scipy
from scipy.integrate import quad
from matplotlib import pyplot as plt

#Trabajare en mks primero
#La idea es hacer una función que me de la matriz sigma (En realidad, me basta con 
#tener un vector de dos elementos pues simga_xx=sigma_yy y sigma_xy=-sigma_yx)
#Las variables son  ...

# constantes
h_barra = 1.054571800*10**-34 #J*s
k_b = 1.3806488*10**-23 #J/k #8.617332478*10**(-5) #eV/K
e = -1.6*10**-19 #C

# parametros fijos
T = 300 #K
mu_c = 0.5 * 1.602177*10**-19 #J
delta = 0
Gamma = 0.11*(10**-3)*(1.602177*10**-19)/(h_barra) #j/hbarra
v_f = 1.02*10**6 #m/s 

# Defino las funciones asi por pasos pq si no redondea mucho

def nF(omega):
    anf=omega-mu_c
    bnf=anf/(k_b*T)
    cnf=np.exp(bnf)
    dnf=cnf+1
    nF=1/dnf
    return(nF)

def M(n):
    aM=delta**2
    bM=2*n*abs(e*B)*h_barra
    cM=v_f**2
    dM= cM*bM+aM
    Mn= np.sqrt(dM) 
    return(Mn)
    

def in_int(g):
    aint = nF(-g)+nF(g)
    bint= omega+2j*Gamma
    bint1 = bint*h_barra
    bint2=bint1**2
    cint= 4* g**2
    dint=cint/bint2
    dint2=1-dint
    eint=aint/dint2
    return(eint)
    

def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])


#%% Sigma xx y sigma yy para B dsitinto de 0

# parametros variables
B = 5 #T = kg/C*s    
    
sigma_xx = []                           #vector con sigma_xx asociados a cada frecuencia
V = np.linspace(3*10**12, 5*10**12, 110)       #frecuencias

for v in V:
    omega = 2*np.pi*v   #2piHz                              #en la formula necesito la frecuencia angular
    n=0                                                    #donde empieza
    omega_fijo = [0, 100]                                  #Estos son para arrancar con el while. Despues no me van a importar. Solo necesito que cumplan la condicion del while y que el segundo elemento no sea cercano a algun valor posible del grafico
    omega_fijo_sinsum = []                                  #vector intermedio en la cuenta
    while abs(omega_fijo[n] - omega_fijo[n+1]) > 1*10**-7.5:         #Fijo la condicion de convergencia
        a= M(n+1)-M(n)
        b= M(n+1)+M(n)
        c= nF(M(n))-nF(M(n+1))+nF(-M(n+1))-nF(-M(n))
        d= nF(-M(n))- nF(M(n+1))+nF(-M(n+1))-nF(M(n))
        e_bis=  h_barra*(omega+2j*Gamma)**2
        a1= 1/a
        a2= a**2/h_barra
        f=a2-e_bis
        g= c/f
        h=a1*g
        b1= 1/b
        b2=b**2/h_barra
        i=b2-e_bis
        j=d/i
        k=b1*j #si pongo b2*j da bien
        in_sum= h+k #hasta aca lo que va dentro de la suma
        a=v_f
        a2=a**2
        b=abs(e*B)*(omega + 2j*Gamma)*h_barra
        c=a2*b
        d=np.pi*1j
        e_bis=c/d
        sinsum = e_bis*in_sum #La expresion sin la sumatoria
        omega_fijo_sinsum.append(sinsum) #Añado la serie hasta n a un vector donde el primer elemento corresponde a n=1, y el segundo a n=2, y asi.
        consum = np.sum(omega_fijo_sinsum)
#        if n % 100 ==0:
#            print(consum)
        omega_fijo.append(consum)
        n=n+1
    print(n, omega)
    sigma_xx.append(omega_fijo[len(omega_fijo)-1])


xx_real = np.real(sigma_xx)
xx_img = np.imag(sigma_xx)
# =============================================================================
# 

# 
# =============================================================================

#%%

sigma_xy = []                           #vector con sigma_xy asociados a cada frecuencia
V = np.linspace(3*10**12, 5*10**12, 110)       #frecuencias

for v in V:
    omega = 2*np.pi*v   #2piHz                              #en la formula necesito la frecuencia angular
    n=0                                                    #donde empieza
    omega_fijo = [0, 100]                                  #Estos son para arrancar con el while. Despues no me van a importar. Solo necesito que cumplan la condicion del while y que el segundo elemento no sea cercano a algun valor posible del grafico
    omega_fijo_sinsum = []                                  #vector intermedio en la cuenta
    while abs(omega_fijo[n] - omega_fijo[n+1]) > 1*10**-8:         #Fijo la condicion de convergencia
        a= M(n+1)-M(n)
        a2=a**2/h_barra
        b= M(n+1)+M(n)
        b2= b**2/h_barra
        c = nF(M(n))-nF(M(n+1))+nF(-M(n+1))-nF(-M(n))
        d2 = h_barra*(omega+2j*Gamma)**2
        e_bis= a2-d2
        e_bis1=1/e_bis
        f_bis = b2-d2
        f_bis1=1/f_bis
        g = e_bis1+f_bis1
        in_sum= c*g #hasta aca lo que esta dentro de la suma
        a=v_f
        a2=a**2
        b= abs(B*e)
        c=a2*b
        d= -c/np.pi
        sinsum = d* in_sum #La expresion sin la sumatoria
        omega_fijo_sinsum.append(sinsum) #Añado la serie hasta n a un vector donde el primer elemento corresponde a n=1, y el segundo a n=2, y asi.
        consum = np.sum(omega_fijo_sinsum)
        omega_fijo.append(consum)
#       if n %1000==0:
#            print(consum)
        n=n+1
    print(n, omega)
    sigma_xy.append(omega_fijo[len(omega_fijo)-1])

xy_real = np.real(sigma_xy)
xy_img = np.imag(sigma_xy)

plt.figure()
plt.plot(V, xy_real)

plt.figure()
plt.plot(V, xy_img)


#%%
np.savetxt('sigma_xx 5T adim ultim.txt', (V, np.real(sigma_xx), np.imag(sigma_xx)))
np.savetxt('sigma_xy 5T adim ultim.txt', (V, np.real(sigma_xy), np.imag(sigma_xy)))


#%% Sigma xx para B igual 0    

sigma_xx0 = []
V = np.linspace(3*10**12, 5*10**12, 110)       #frecuencias


for v in V:
    omega = 2*np.pi*v   #2piHz                              #en la formula necesito la frecuencia angular
    n=0   
    limite=0                                                 #donde empieza
    omega_fijo = [0, 100]                                  #Estos son para arrancar con el while. Despues no me van a importar. Solo necesito que cumplan la condicion del while y que el segundo elemento no sea cercano a algun valor posible del grafico                               #vector intermedio en la cuenta
    while abs(np.real(omega_fijo[n] - omega_fijo[n+1])) > 1*10**-14:         #Fijo la condicion de convergencia
        a= 1 #e**2
        b= 4 #/h_barra
        c=a/b
        d=4j*c
        e_bis=d/np.pi
        f= omega+2j*Gamma
        f1=1/f
        g= 1/k_b*T
        g1=2/g
        h = -mu_c*g1
        h1=np.exp(h)
        h2=np.log(1+h1)
        i=h2*g1
        i1=mu_c+i
        k = complex_quadrature(in_int, 0, limite)[0]
        m = k+i1
        m1= m/(h_barra*f)
        omega_fijo_n=e_bis*m1 #hasta aca la formula
        omega_fijo.append(omega_fijo_n)
        limite=(n+1)*100
        n=n+1
    sigma_xx0.append(omega_fijo[len(omega_fijo)-1])

np.savetxt('sigma_xx 0T adim desde3.txt', (V, np.real(sigma_xx0), np.imag(sigma_xx0)))

#%%#Grafico 
xx_real0 = np.real(sigma_xx0)
xx_img0 = np.imag(sigma_xx0)

plt.figure()
plt.plot(V, xx_real0)

plt.figure()
plt.plot(V, xx_img0)

omegon = 2*np.pi*V*h_barra/mu_c

plt.figure()
plt.plot(omegon, xx_real0)

plt.figure()
plt.plot(omegon, xx_img0)





