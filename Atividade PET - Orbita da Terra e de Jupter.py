from numpy import array
from pylab import plot, show
from vpython import sphere,rate,vector,color, canvas

# Dados tirados do site: https://radiojove.gsfc.nasa.gov/education/educ/jupiter/basics/jfacts.htm

def f(r,s,t):    #EDOs que regem o movimento do planeta r influenciado gravitacionalmente pelo sol e pelo planeta s.
    sx,svx,sy,svy,rx,rvx,ry,rvy  = s[0],s[1],s[2],s[3],r[0],r[1],r[2],r[3]
    R = (rx**2 + ry**2)**(1/2)     #Distância entre o planeta r e o sol.
    RS = ((rx-sx)**2 + (ry-sy)**2)**(1/2)     #Distância entre o planeta r e o planeta s.
    
    #Equações diferenciais de primeira ordem
    drx = rvx     #drx/dt
    dry = rvy     #dry/dt
    drvx = -Gm*(rx-sx)/(RS**3) -GM*rx/(R**3)      #drvx/dt
    drvy = -Gm*(ry-sy)/(RS**3) -GM*ry/(R**3)      #drvx/dt
    
    return array([drx,drvx,dry,drvy])

def RK4(r,s,t):     #Método Runge Kutta 4ª ordem
    k1 = h*f(r,s,t)
    k2 = h*f(r+k1/2,s,t+h/2)
    k3 = h*f(r+k2/2,s,t+h/2)
    k4 = h*f(r+k3,s,t+h)
    return r + (1/6)*(k1+2*k2+2*k3+k4)

# Parâmetros do cálculo:
a = 0.0     #Tempo inicial.
b = 12*3.1536*10**7     #Tempo final (12x Periodo da orbita da Terra). 
N = 5*10**4     #Numero de passo da simulação.
h = (b-a)/N     #Intervalo de tempo de cada passo da simulação.
t = a     #Parânemtro de calculo que representa o tempo(s)
n = 0     #Parâmetro auxilixar
G = 6.6738*10**(-11)     #Constante da gravitação universal
M = 1.9891*10**30     #Massa do sol(kg)
m1 = 5.9736*10**24     #Massa da Terra(kg)
m2 = 1.8986*10**27     #Massa de Jupiter(kg)
GM = G*M
Gm1 = G*m1
Gm2 = G*m2


temp = []     #Tempo (t minúsculo).
Tx = []     #Componetes da posição e velocidade da Terra nas coordenadas x e y em todos os instantes de tempo.
Ty = []
Tvx = []
Tvy = []
Jx = []     #Componetes da posição e velocidade de Jupiter nas coordenadas x e y em todos os instantes de tempo.
Jy = []
Jvx = []
Jvy = []

#OBS: Comprimentos em metro e velocidades em metro por segundo.

Tx0,Tvx0,Ty0,Tvy0 = 1.4709*10**11,0.0,0.0,3.029*10**4     #Posição e velocidade inicial da Terra (no perielio).
Jx0,Jvx0,Jy0,Jvy0 = 7.4052*10**11,0.0,0.0,1.372*10**4      #Posição e velocidade inicial de Jupiter (no perielio).
T = array([Tx0,Tvx0,Ty0,Tvy0])     #Posição e velocidade inicial da Terra em um array (T maiúsculo).
J = array([Jx0,Jvx0,Jy0,Jvy0])     #Posição e velocidade inicial de Jupiter em um array (J maiúsculo).

#OBS: Para facilitar os cálculos vamos considerar que tanto Jupiter quanto a Terra partem dos seus respectivos 
#   perielios no instante inicial. Alem disso a posição desses pontos são relativimente próximas, o que pode não
#   condizer com a realidade.

# Resolvendo EDO:
while(t<=b):
    temp.append(t)     #Salvando o tempo atual.
    Tx.append(T[0])     #Salvando a posição e velocidade da Terra no tempo t atual.
    Tvx.append(T[1])
    Ty.append(T[2])
    Tvy.append(T[3])
    Jx.append(J[0])     #Salvando a posição e velocidade de Jupiter no tempo t atual.
    Jvx.append(J[1])
    Jy.append(J[2])
    Jvy.append(J[3])
    
    #Calculando a posição e velocidade da Terra no próximo passo.
    Gm = Gm2
    S = array([Jx[n],Jvx[n],Jy[n],Jvy[n]])
    T = RK4(T,S,t)
    
    #Calculando a posição e velocidade de Jupiter no próximo passo.
    Gm = Gm1
    S = array([Tx[n],Tvx[n],Ty[n],Tvy[n]])
    J = RK4(J,S,t)
    
    n+=1
    t+=h     #Evoluindo o tempo t atual
    
#plot(xT,yT)     #Plotando a orbita da Terra
#plot(xJ,yJ)     #Plotando a orbita de Jupiter
#show()

scene = canvas(width = 900, height = 500, center = vector(0,0,0) ,background=color.black)

#T1 = sphere(pos=vector(0,0,0),radius=e*6.3781*10**6,color = color.blue)     #Terra
#J1 = sphere(pos=vector(0,0,0),radius=e*7.1492*10**7,color = color.orange)     #Jupiter
#S1 = sphere(pos=vector(0,0,0),radius=e*6.96340*10**8,color = color.white)     #Sol

T1 = sphere(pos=vector(0,0,0),radius=10**10,color = color.blue)     #Terra
J1 = sphere(pos=vector(0,0,0),radius=10**10,color = color.orange)     #Jupiter
S1 = sphere(pos=vector(0,0,0),radius=10**10,color = color.white)     #Sol

for i in range(n):
    rate(700)
    T1.pos = vector(Tx[i],Ty[i],0)
    J1.pos = vector(Jx[i],Jy[i],0)