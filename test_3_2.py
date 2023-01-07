import math
import cmath
import numpy as np
import sympy as sp
import os

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

def DJ1J0 (Z):
    J0 = sp.besselj(0,Z)
    J1 = sp.besselj(1,Z)

    return J1/J0
filef = open('./G.txt', 'a', encoding='UTF-8')
bmax = 100

PI = math.pi                    # 正直、pythonだとこの定義は必要ないかも？
I = complex(0.0,1.0)            # 上に同じ

ME = Ms = Mf = Mff = Msta = np.zeros((2,2),dtype=complex)


# 実験条件
R = 0.34 * (10**(-3))           # 蓄熱器の流路半径R粘性
Pm = 101.0 * (10**3)            #平均気圧
Tc=293                          #蓄熱器低温端の温度
TH=293*(2.3)                    #蓄熱器高温端の温度
Lstack=20* (10**(-3))
f=100.0                         #周波数

# 計算の入口--------------------------------------------------
w = 2.0 * f * PI

Tm = Tc

rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)

Pin=complex(1.08371789371,-0.19130178936)     
Uin=complex(1.777853001082 * (10**(-3)),8.733481219974 * (10**(-5)))

Iin = 0.5 * (Pin.conjugate() * Uin).real

Mff = ME = np.identity(2,dtype=complex) 
# 計算の入口終わり----------------------------------------------

DelX=Lstack/bmax        #管の刻み幅
dTm=(TH-Tc)/Lstack      #温度勾配

for b in range(1,bmax):
    Tm = Tc+dTm*DelX*b
    rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)
    wtu =((R**2)*w)/(2.0*Nyu ) #このあたりの記号の使い方はは富永著の熱音響工学の基礎を参照
    wta =((R**2)*w)/(2.0*Alpha )

    Ya = complex(1.0,1.0)*complex(cmath.sqrt(wta),0.0) # ここ*complex()って必要？
    Yu = complex(1.0,1.0)*complex(cmath.sqrt(wtu),0.0)

    Xa = ((2.0/I)/Ya)*DJ1J0(I*Ya)
    Xu = ((2.0/I)/Yu)*DJ1J0(I*Yu)

    Ms[0,0] = 0.0
    Ms[0,1] = (-I*w*rhom)/(1.0-Xu)
    Ms[1,0] = -I*w*(1.0+(Gam-1.0)*Xa)/(Gam*Pm)
    Ms[1,1] = (Xa-Xu)*dTm/(1.0-Xu)/(1.0-Pr)/Tm

    Mf[0,0] = ME[0,0] + DelX*Ms[0,0]
    Mf[0,1] = ME[0,1] + DelX*Ms[0,1]
    Mf[1,0] = ME[1,0] + DelX*Ms[1,0]
    Mf[1,1] = ME[1,1] + DelX*Ms[1,1]

    Msta = np.matmul(Mf,Mff)    # Msta = Mf @ Mff

    Mff = Msta

    # numpyでは数値とリストの積'*'で、定数倍になる。
    # 普通のリストだと定数倍は定義されていないので注意

    P = Pin*Mff[0,0] + Uin*Mff[0,1]
    U = Pin*Mff[1,0] + Uin*Mff[1,1]
    In = 0.5 * (P*(U.conjugate())).real
    print( abs(P), abs(U), (In/Iin), file=filef)

Pout = Mff[0,0]*Pin + Mff[0,1]*Uin
Uout = Mff[1,0]*Pin + Mff[1,1]*Uin
Iout = 0.5* ((Pout.conjugate() * Uout).real)

print('増幅率', Iout/Iin, '周波数',f)
print('出口の音響インピーダンス',Pout/Uout)
print('出口での特性音響インピーダンス',rhom*Ss)

