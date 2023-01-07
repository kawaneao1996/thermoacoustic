# Sigma_A_kappaを求めるところが未完成
# Haを求めるところが未実装
import math
import cmath
import numpy as np
# import sympy as sp
from scipy.special import jv
import os
# mmax=1    #エンタルピーの値を探索する用
bmax=100  #スタック内を刻む用

Pi = math.pi                    # 正直、pythonだとこの定義は必要ないかも？
I = complex(0.0,1.0)            # 上に同じ
ME = np.identity(2,dtype=complex)
Ms = Mf = np.zeros((2,2),dtype=complex)
    
R = 0.34 * (10**(-3))           # 蓄熱器の流路半径R
Pm = 101.0 * (10**3)            #平均気圧
    
Ta=18 + 273                          #蓄熱器一端の温度
Tb=278 + 273                    #蓄熱器他端の温度
Lstack=35 * (10**(-3))
f= 41                        #周波数
w = 2.0 * f * Pi
DelX=Lstack/bmax        #管の刻み幅    
Sigma_A_kappa = 0 # ここに熱伝導の項を入れる。DelTm/DelXを計算する分母の計算に使う。


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
    # なぜかこのベッセル関数がsympyで計算できなくなった
    # J0 = sp.besselj(0,Z)
    # J1 = sp.besselj(1,Z)
    J0 = jv(0, Z)
    J1 = jv(1, Z)
    # print("Z" , Z)
    #print("J0" , J0)
    #print("J1" , J1)
    rslt = J1 / J0
    if (J0 == complex(float("nan"),float("nan"))):
        return 0
    # print("J1/J0",J1/J0)
    return rslt


def ret_Tm (Ha,Ta,Pin,Uin,flag=False):
    global I,w,bmax,R,Pm,Lstack,f # 条件値で書き直されないもの
    H = Ha
    Tm = Ta
    P = Pin
    U = Uin

    for b in range(1,bmax):
        ##########################
        # 蓄熱器内の計算スタート #
        ##########################
        
        rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)
        
        wtu =( R ** 2) * w / (2.0 * Nyu)
        
        wta =  (R ** 2) * w / (2.0 * Alpha)
        
        Ya = complex(1.0,1.0) * complex(cmath.sqrt(wta), 0.0)
        Yu = complex(1.0,1.0) * complex(cmath.sqrt(wtu), 0.0)

        YaI = I * Ya
        YuI = I * Yu
        # print("YaI",YaI)
        # print("DJ1J0",DJ1J0(YaI))
        Xa = 2.0 / I / Ya * DJ1J0(YaI)
        Xu = 2.0 / I / Ya * DJ1J0(YuI)

        if (Xa == 0 and Xu == 0):
            return "zero_devide"
        # 熱伝導の項を入れるなら、流路半径Rが必要（上田先生のプログラムでは分子分母をA1＝Rで割って省略されている）
    # pythonでは改行は'\'で行う

        C1 = R * rhom * Cp * (abs(U)**2)
        C2 = 2.0 * w * (1.0 - Pr * Pr) * abs(1.0 - Xu) ** 2
        # print(Xa)
        # print(Pr)
        # print(Xu)
        C3 = (Xa + Pr * Xu.conjugate())
        # print(C3)
        C3 = C3.imag
        C = C1 / C2 * C3 - Sigma_A_kappa

        # 熱伝導の項をSigma_A_kappaとした
        
        DelTm = (H - 0.5 * (P * U.conjugate() * (1.0 - (Xa - Xu.conjugate()) / ((1 + Pr) * (1 - Xu.conjugate())) ) ).real) / C
        
        Ms[0][0] = 0.0
        Ms[0][1] = (- I * w * rhom)/ (1.0 - Xu)
        Ms[1][0] = - I * w * (1.0 + (Gam - 1.0) * Xa)/(Gam * Pm)
        Ms[1][1] = (Xa - Xu) * DelTm / (1.0 - Xu) / (1.0 - Pr) / Tm # /DelTm -> /Tm
        
        Mf = ME + DelX * Ms # pythonではnp.arrayに対して'*'で成分の定数倍、'+'で行列の成分同士の足し算ができる
        
        
        P_tmp = P
        U_tmp = U
        
        P = P_tmp * Mf[0,0] + U_tmp * Mf[0,1]
        U = P_tmp * Mf[1,0] + U_tmp * Mf[1,1]
        
        Tm = Tm + DelTm * DelX
        In = 0.5 * (P * U.conjugate()).real
        if flag:
            print(DelX*b,abs(P),abs(U),file=fileA)
            print(DelX*b,Tm,file=fileT)
            print(DelX*b,H,In,H-In,file=filee)
        
    return Tm,P,U

def judge_Ha(Ha,Ta,Tb,Pin,Uin,threshold):
    T_tmp,_,_ = ret_Tm(Ha,Ta,Pin,Uin)
    diff = abs(Tb - T_tmp)
    if (T_tmp == "zero_devide"):
        return False
    if (diff < threshold):
        return True
    else:
        return False

def main():

    I = complex(0.0,1.0)
    fileA = open('./Acoustic_filed.txt', 'a', encoding='UTF-8') # 位置、圧力振幅，仕事流
    fileT = open('./Temperature.txt', 'a', encoding='UTF-8')    # 位置，時間平均温度
    filee = open('./enegy_flow.txt', 'a', encoding='UTF-8') # 位置，エンタルピー流，仕事流，熱流
    
    
    ##############
    # 入口の計算 #
    ##############
    

    Tm = Ta
    
    rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)
    
    ####################
    # 入口の圧力と流速 #
    ####################
    
    Pin = complex(3.4 * (10**3), 0.0)
    Uin = complex(3.4 * (10**3)/(9.0 * rhom * Ss), 0.0)
    
    Iin = 0.5 * (Pin.conjugate() * Uin).real
    
    ######################
    # エンタルピーの計算 #
    ######################
    Ha = 0.12 * Iin

    Mff = ME
    DelX=Lstack/bmax        #管の刻み幅
    # ここでHaは求まっている
    P = Pin
    U = Uin
    H = Ha
    for i in range(1,bmax):
        rhom,Mu,Kappa,Gam,Cp,Nyu,Alpha,Pr,Ss = bussei(Tm,Pm)
        wtu =(R**2)* w /(2.0 * Nyu)             
        wta =(R**2)* w /(2.0 * Alpha )
      
        Ya =complex(1,1)*complex(cmath.sqrt(wta),0)
        Yu =complex(1,1)*complex(cmath.sqrt(wtu),0)            
      
        Xa =2.0/I/Ya*DJ1J0(I*Ya)
        Xu =2.0/I/Yu*DJ1J0(I*Yu)
        C1 = R * rhom * Cp * (abs(U)**2)
        C2 = 2.0 * w * (1.0 - Pr * Pr) * abs(1.0 - Xu) ** 2
        # print(Xa)
        # print(Pr)
        # print(Xu)
        C3 = (Xa + Pr * Xu.conjugate())
        # print(C3)
        C3 = C3.imag
        C = C1 / C2 * C3 - Sigma_A_kappa

        # 熱伝導の項をSigma_A_kappaとした
        
        DelTm = (H - 0.5 * (P * U.conjugate() * (1.0 - (Xa - Xu.conjugate()) / ((1 + Pr) * (1 - Xu.conjugate())) ) ).real) / C
        
        Ms[0][0] = 0.0
        Ms[0][1] = (- I * w * rhom)/ (1.0 - Xu)
        Ms[1][0] = - I * w * (1.0 + (Gam - 1.0) * Xa)/(Gam * Pm)
        Ms[1][1] = (Xa - Xu) * DelTm / (1.0 - Xu) / (1.0 - Pr) / Tm # /DelTm -> /Tm
        
        Mf = ME + DelX * Ms # pythonではnp.arrayに対して'*'で成分の定数倍、'+'で行列の成分同士の足し算ができる
        
        
        P_tmp = P
        U_tmp = U
        
        P = P_tmp * Mf[0,0] + U_tmp * Mf[0,1]
        U = P_tmp * Mf[1,0] + U_tmp * Mf[1,1]
        
        Tm = Tm + DelTm * DelX
        In = 0.5 * (P * U.conjugate()).real

        Pout = P
        Uout = U
        Iout = 0.5 * ((Pout.conjugate() * Uout).real)

    print("想定温度", Tb)
    print("計算の結果得られた出口温度",Tm)
    print("差",Tb - Tm)
    print("増幅量", In - Iin)
    print("出口における熱流", H - In)
    print("効率", -(In - Iin)/(H - In))


if __name__ == '__main__':
    main()



