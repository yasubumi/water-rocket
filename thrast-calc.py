
# coding: utf_8
# water rocket thrust calcuration
# Date: 2018-05-03
# Version: 0.3
'''
   Copyright 2017 Yasubumi KANAZAWA (camelinsect@wings2fly.jp)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
'''

import math
import sys

#----- 設定 -----
St = (9.5/2)**2*math.pi/1000.0**2 #ノズル出口断面積[m^2]
dt = 0.001 #時間ステップ[s]
m=0.134 # 機体質量 [kg]
P0 = 686000.0 # 初期タンク内圧[Pa] (7kgf/cm^2 = 100 PSI)
V0 = 0.670/1000 # 初期タンク内水体積(0.5L) [m^3]
T0 = 300.0 #初期タンク内空気温度[K]
dLpet = 200.0 #タンク長さ延長[mm]:
Cd = 0.75 # 機体の抗力係数
A = (0.09/4)**2*math.pi #機体の断面積[m^2]

#----- 定数 -----
g=9.8 #重力加速度[m/s^2]
PA = 101300.0 #大気圧[Pa]
gamma = 1.4 #比熱比
rho = 1000.0 #水密度 [kg/m^3]
rho_air = 1.293 #空気密度[kg/m^3]

#----- 変数 -----
X0 = 0.0 #初期水面位置[m]
S0 = 0.0 #タンク断面積[m^2]
Vair0 = 0.0 #初期タンク内空気体積[m^3]
t = 0.0 #時間
P = 0.0 #タンク内圧[Pa]
V = 0.0 #タンク内水体積[m^3]
S = 0.0 #水面位置でのタンク断面積[m^2]
Vair = 0.0 #タンク内空気体積[m^3]
T = 0.0 #タンク内空気温度[K]
v = 0.0 #水噴出速度[m/s]
x = 0.0 #水面位置[m]
F = 0.0 #推力[N]
M = 0.0 #質量[kg]
u = 0.0 #ロケット速度[m/s]
Ft = 0.0 #力積[Nsec]
H_tmp = 0.0 #概算到達高度[m]
z = 0.0 #ロケット高度[m]


#----- ペットボトル関連 -----
Xpet = [65.0,245.0+dLpet,275.0+dLpet,301.0+dLpet] # ペットボトル形状が変化する位置(ボトル底面からの距離)
Vp = [0,0,0] # ペットボトルの体積
Vp[0] = 0.0094393/1000        # Xpet[3]まで水を入れたとき…出口部
Vp[1] = Vp[0]+1.0/3.0*math.pi*(0.01075**2+0.01075*0.088/2+(0.088/2)**2)*0.030  # Xpet[2]まで水を入れたとき…徐変部
Vp[2] = Vp[1]+(0.088/2)**2*math.pi*(Xpet[1]-Xpet[0])/1000  # Xpet[2]まで水を入れたとき…円筒部
Vpet = Vp[2]+0.3602150/1000   # Xpet[0]まで水を入れたとき…ペットボトル全体の体積

apet=[0,0,0,0]
apet[0]=((0.044-0.01075)/0.030)**2
apet[1]=3*0.01075*apet[0]
apet[2]=3*0.01075**2
apet[3]=-3/math.pi

#ペットボトル内の水量から水面位置を求める
def calc_x(V):
    dV=0.0
    dx=0.0
    x=0.0
    if V>Vpet:
        print "ペットボトル容量オーバー"
        return -1
    elif V>Vp[2]:
        print "水量制限オーバー"
        return -1
    elif V>Vp[1]:
        #print "円筒部"
        dV=V-Vp[1]
        dx=dV/(0.044**2*math.pi)
        x=Xpet[1]/1000-dx
    elif V>Vp[0]:
        #print "拡大部"
        dV=V-Vp[0]
        dx=0.001
        while True:
            dx1 = dx -(apet[0]*dx**3+apet[1]*dx**2+apet[2]*dx+apet[3]*dV)/(3*apet[0]*dx**2+2*apet[1]*dx+apet[2])
            if abs(dx1-dx) < 0.0000001: 
                dx=dx1 
                break 
            dx=dx1 
        x=Xpet[2]/1000-dx 
    elif V>0:
        #print "出口部"
        dV=V
        dx=dV/(0.01075**2*math.pi)
        x=Xpet[3]/1000-dx
    else:
        print "水体積がマイナス"
        return -1
    return x

#水面位置から水面の面積を求める
def calc_S(x):
    S=0.0
    r=0.0
    dx=0.0
    if x>Xpet[3]/1000:
        print "水面位置がボトルの外"
        return -1
    elif x<0.0:
        print "水面位置がマイナス"
        return -1
    elif x<Xpet[0]/1000:
        print "水容量オーバー"
        return -1
    elif x<Xpet[1]/1000:
        #print "円筒部"
        r=0.044
    elif x<Xpet[2]/1000:
        #print "徐変部"
        dx=x-Xpet[2]/1000
        r=0.01075+(0.044-0.01075)/0.030*dx
    elif x>=Xpet[2]/1000:
        #print "出口部"
        r=0.01075
    S=r**2*math.pi
    return S

#----debug----
for i in range(3):
    print "Vp["+str(i)+"]="+str(Vp[i])
print "Vpet"+str(Vpet)
for i in range(4):
    print "Xpet["+str(i)+"]="+str(Xpet[i])

#----end debug----


#----- 初期化 -----
x0 = calc_x(V0) #初期水面位置[m]
print "x0=" + str(x0) 
if x0 == -1:
    print "初期水面位置エラー"
    sys.exit(-1)
S0 = calc_S(x0) # タンク断面積[m^2]
if S0 == -1:
    print "初期水面面積エラー"
    sys.exit(-1)
Vair0 = Vpet-V0 # 初期タンク内空気体積[m^3]
t=0.0 # 時間
P=P0 # タンク内圧
V=V0 # タンク内水体積[m^3]
S=S0 # 水面位置のタンク断面積
Vair = Vair0 # タンク内空気体積[m^3]
T=T0 # タンク内空気温度[K]
v=0.0 # 噴出速度
x=x0 # 水面位置
F=0.0 # 推力
M=m+V0*rho # 全備質量

#----- 関数 -----
# Vの変数分離型 微分方程式
def f_V(V):
    x=calc_x(V)
    S=calc_S(x)
    dVdt = St* math.sqrt(2/rho * (P0*(Vair0/(Vpet-V))**gamma -PA))
    return dVdt

# ルンゲクッタでVを求める
def calc_Vu(V,u):
    k0=[0.0,0.0,0.0,0.0]
    k1=[0.0,0.0,0.0,0.0]

    k0[0] = f_V(V)
    k1[0] = f_u(V,u)

    k0[1] = f_V(V+dt*0.5*k0[0])
    k1[1] = f_u(V+dt*0.5*k0[0], u+dt*0.5*k1[0])

    k0[2] = f_V(V+dt*0.5*k0[1])
    k1[2] = f_u(V+dt*0.5*k0[1], u+dt*0.5*k1[1])

    k0[3] = f_V(V+dt*k0[2])
    k1[3] = f_u(V+dt*k0[2], u+dt*0.5*k1[2])

    V1 = V-dt/6.0 *(k0[0]+2.0*k0[1]+2.0*k0[2]+k0[3])
    u1 = u+dt/6.0 *(k1[0]+2.0*k1[1]+2.0*k1[2]+k1[3])
    return V1,u1

# ロケット速度uの微分方程式
def f_u(V,u):
    dudt=f_F(V)/f_M(V)-g-(0.5*rho_air*Cd*A)/f_M(V)*u**2
    return dudt

# 圧力Pの式
def f_P(V):
    ans = 0.0
    ans = P0*(Vair0/(Vpet-V))**gamma
    return ans

# 水の噴出速度vの式
def f_v(V):
    ans=0.0
    ans=(2/rho*(f_P(V)-PA))**0.5
    return ans

# 推進力Fの式
def f_F(V):
    ans = 0.0
    ans = rho*St*f_v(V)**2
    return ans

# 質量Mの式
def f_M(V):
    ans = 0.0
    ans = m+rho*V
    return ans


#----- メインプログラム -----
# 初期値確認
print "初期設定値"
print "時間ステップ = "+str(dt)+"[sec]"
print "ノズル断面積 = "+str(St*10**6)+"[mm^2]"
print "初期水量 = "+str(V0*1000)+"[L]"
print "機体質量 = "+str(m)+"[kg]"
print "初期タンク内圧 = "+str(P/100)+"[hPa]"
print "初期タンク内温度 = "+str(T0-273.15)+"[℃]"
print "\n"

# 水噴射ステージ
print "時間(t)[sec] 水面位置(x)[m] 圧力(P)[Pa] 噴出速度(v)[m/s] 水面面積(S)[m^2] ロケット質量(M)[kg] 推力(F)[N] 水体積(V)[L] 力積(累計)[Nsec] ロケット速度(u)[m/s] ロケット位置(z)[m]"
#初期(t=0)
v=f_v(V)
F=f_F(V)
Ft=F*dt
print str(t),
print str(x),
print str(P),
print str(v),
print str(S),
print str(M),
print str(F),
print str(V),
print str(Ft),
print str(u),
print str(z)
t+=dt
while True:
    V1,u1=calc_Vu(V,u)
    if V1<0: # ここで残りの水噴射の計算を実施 if V>Vp[0]:
            # 水面が出口部に達していなければエラーで終了
        if V>Vp[0]:
            print "時間当たりの水体積変化が大きすぎます。時間の刻みを小さくしてください"
            break
        else:
            # 残りは噴出速度vを一定として排出時間とタンク内圧力変化を計算。
            x=Xpet[3]/1000
            P=f_P(0)
            if P-PA<0.1:
                print "タンク内圧が大気圧まで下がりました。初期圧力を上げてください"
                break
            S=calc_S(x)
            M=m
            v=f_v(0)
            h=V/(v*St)
            t+=h
            F=f_F(0)
            V=0
            Ft=Ft+F*dt
            dt0 = dt
            dt=h
            u1 = u+(F/M-g-0.5/M*rho_air*Cd*A*u**2)*dt
            z=z+(u+u1)*0.5*dt
            u=u1
            dt=dt0
            print str(t),
            print str(x),
            print str(P),
            print str(v),
            print str(S),
            print str(M),
            print str(F),
            print str(V),
            print str(Ft),
            print str(u),
            print str(z)
            print "水噴出ステージ正常終了"
            break
    z=z+(u+u1)*0.5*dt
    V=V1
    u=u1
    x=calc_x(V)
    Vair=Vpet-V
    P=f_P(V)
    if P-PA<0.1:
        print "タンク内圧が途中で大気圧まで下がりました。初期圧を上げてください"
        break
    v=f_v(V)
    S=calc_S(x)
    M=f_M(V)
    F=f_F(V)
    Ft=Ft+F*dt
    print str(t),
    print str(x),
    print str(P),
    print str(v),
    print str(S),
    print str(M),
    print str(F),
    print str(V),
    print str(Ft),
    print str(u),
    print str(z)
    t+=dt

