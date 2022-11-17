import math
import cmath
import numpy as np
import sympy as sp
import os

flag_forward = False #進行波が正の方向に進行する/False だと負の方向に進む

# program tubeの内容
##############
# 変数の宣言 #
##############
pi = math.pi
I = complex(0.0,1.0)
N : int  = 0
Nmax : int = 100 #管の分割数

# 管の形状など
Pm = 101.0 * 10**3 #mean pressure
r = 10.5 * 10**(-3) #radious of tube
L = 1.00    #length of tube
A = r * r * pi #cross-sectional area of tube
TH = 300.0  #mean temperature
f =  148.5  #frequency
omega = 2.0 * pi * f   #angular frequency
if flag_forward:
    dx = - L / float(Nmax)
else:
    dx =  L / float(Nmax)

############
# 境界条件 #
############

#P0:pressure at tube-end 閉端の場合は有限値、開端の場合は０
P0 = complex(0.78 * 10**3, 0.0)

#U0 : velocity at tube-end 閉端の場合は０、開端の場合は有限地
U0 = complex(0.0, 0.0) 

Tm = TH

def bussei (Tm,Pm):
    # C密度
    rhom = (4.1461 
        -0.019705*Tm
        +4.7881*(10**(-5))*Tm*Tm
        -6.2842*(10**(-8))*Tm*Tm*Tm
        +4.2339*(10**(-11))*Tm*Tm*Tm*Tm
        -1.1472*(10**(-14))*Tm*Tm*Tm*Tm*Tm)*Pm/(101.0 * 10**(3))
    # C粘性係数
    Mu = -7.0031*(10**(-6))               
    +0.14018*(10**(-6))*Tm           
    -0.00026899*(10**(-6))*Tm*Tm      
    +3.5727*(10**(-13))*Tm*Tm*Tm      
    -2.4676*(10**(-16))*Tm*Tm*Tm*Tm   
    +6.8396*(10**(-20))*Tm*Tm*Tm*Tm*Tm
    # C熱伝導度
    Kappa=4.9135*(10**(-3))
    +0.053336*(10**(-3))*Tm
    +0.00013329*(10**(-3))*Tm*Tm
    -3.4467*(10**(-10))*Tm*Tm*Tm
    +3.4584*(10**(-13))*Tm*Tm*Tm*Tm
    -1.2551*(10**(-16))*Tm*Tm*Tm*Tm*Tm

    # C比熱比
    Gam=1.3803
    +0.00020301*Tm
    -5.1868*(10**(-7))*Tm*Tm
    +2.4927*(10**(-10))*Tm*Tm*Tm
    +1.0045*(10**(-13))*Tm*Tm*Tm*Tm
    -7.8378*(10**(-17))*Tm*Tm*Tm*Tm*Tm
    # C定圧比熱
    Cp=1057.9
    -0.3684*Tm
    +0.00063128*Tm*Tm
    +3.4653*(10**(-7))*Tm*Tm*Tm
    -8.7851*(10**(-10))*Tm*Tm*Tm*Tm 
    +3.6441*(10**(-13))*Tm*Tm*Tm*Tm*Tm 
    # C動粘性係数     
    Nyu =Mu /rhom
    # C熱拡散係数
    Alpha = Kappa / rhom / Cp
    # Cプラントル数
    Pr =Mu / Kappa *Cp
    # C音速
    Ss =cmath.sqrt(Gam * Pm /rhom )

    return rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss

rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)

wta = (r * r) * omega / (2.0 * Alpha)

wtu = (r * r) * omega / (2.0 * Nyu)

Ya = complex(1.0,1.0) * complex(cmath.sqrt(wta), 0.0)
Yu = complex(1.0,1.0) * complex(cmath.sqrt(wtu), 0.0)

def DJ1J0 (Z):
    J0 = sp.besselj(0,Z)
    J1 = sp.besselj(1,Z)

    return J1/J0

Xa = 2.0 / I / Ya * DJ1J0(I * Ya)
Xu = 2.0 / I / Yu * DJ1J0(I * Yu)

K0 = omega/Ss

K = K0 * cmath.sqrt(( 1.0 + (Gam - 1.0) * Xa) / (1.0 - Xu))

Mtube = np.array([[0,0],[0,0]], dtype=complex)
Mtube1 = np.array([[0,0],[0,0]], dtype=complex)

f = open('./myfile.txt', 'a', encoding='UTF-8')


for N in range(1,Nmax):
  X = N * dx

  Mtube[0][0]=cmath.cos(K * X)
  Mtube[0][1]=-(I * omega * rhom)/(K * (1-Xu)) *  cmath.sin(K * X)
  Mtube[1][0]=(K * (1-Xu))/(I * omega * rhom) *  cmath.sin(K * X)  
  Mtube[1][1]=cmath.cos(K * X)

  P = P0 * Mtube[0][0] + U0 * Mtube[0][1]
  U = P0 * Mtube[1][0] + U0 * Mtube[1][1]

  Pamp = abs(P)
  Uamp = abs(U)

  Pphi = math.atan2(P.imag, P.real) # = cmath.phase(P)
  Uphi = math.atan2(U.imag, U.real)

  phi = Uphi - Pphi

  tmp = P * (U.conjugate())
  
  W = 0.5 * tmp.real

  print(X,Pamp,Uamp,Pphi,W,file=f)
  # datalist = [X,Pamp,Uamp,Pphi,W]
  # line = str(datalist)
  # f.writelines(line)
  
f.close()
