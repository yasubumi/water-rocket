# coding: utf_8
# water rocket thrust calcuration
# Date: 2017-05-12
# Version: 0.1
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
St = (7.7/2)**2*math.pi/1000.0**2 #ノズル出口断面積[m^2]
dt = 0.05 #時間ステップ[s]
m=0.134 # 機体質量 [kg]
P0 = 686000.0 # 初期タンク内圧[Pa]
V0 = 1.000/1000 # 初期タンク内水体積(0.5L) [m^3]
T0 = 300.0 #初期タンク内空気温度[K]

#----- 定数 -----
g=9.8 #重力加速度[m/s^2]
PA = 101300.0 #大気圧[Pa]
gamma = 1.4 #比熱比
rho = 1000.0 #水密度 [kg/m^3]

#----- ペットボトル関連 -----
Xpet = [65.0,187.0,200.0,240.0,280.0,306.0] # ペットボトル形状が変化する位置(ボトル底面からの距離)
Vp = [0,0,0,0,0] # ペットボトルの体積
Vp[0] = 0.0094393/1000        # Xpet[5]まで水を入れたとき…出口部
Vp[1] = Vp[0]+0.1379391/1000  # Xpet[4]まで水を入れたとき…2次関数部
Vp[2] = Vp[1]+0.2546366/1000  # Xpet[3]まで水を入れたとき…拡大部
Vp[3] = Vp[2]+0.0864325/1000  # Xpet[2]まで水を入れたとき…縮小部
Vp[4] = Vp[3]+0.7761305/1000  # Xpet[1]まで水を入れたとき…円筒部（ここまでを水を入れる許容量とする）
Vpet = Vp[4]+0.3602150/1000   # Xpet[0]まで水を入れたとき…ペットボトル全体の体積

#ペットボトル内の水量から水面位置を求める
def calc_x(V):
    dV=0.0
    dx=0.0
    x=0.0
    if V>Vpet:
        print "ペットボトル容量オーバー"
        return -1
    elif V>Vp[4]:
        print "水量制限オーバー"
        return -1
    elif V>Vp[3]:
        #print "円筒部"
        dV=V-Vp[3]
        dx=dV/(0.045**2*math.pi)
        x=0.187-dx
    elif V>Vp[2]:
        #print "拡大部"
        dV=V-Vp[2]
        dx=0.007
        while True:
            dx1 = dx -(0.024785741*dx**3-0.022716131*dx**2+0.006939778*dx-dV)/(0.074357223*dx**2-0.045432263*dx+0.006939778)
            if abs(dx1-dx) < 0.0000001:
                dx=dx1
                break
            dx=dx1
        x=0.2-dx
    elif V>Vp[1]:
        #print "縮小部"
        dV=V-Vp[1]
        dx=0.01
        while True:
            dx1 = dx -(0.010471976*dx**3+0.013508848*dx**2+0.005808805*dx-dV)/(0.031415927*dx**2+0.027017697*dx+0.005808805)
            if abs(dx1-dx) < 0.0000001:
                dx=dx1
                break
            dx=dx1
        x=0.24-dx
    elif V>Vp[0]:
        #print "2次関数部"
        dV=V-Vp[0]
        dV=0.1379391/1000 - dV
        dx=10
        while True:
            dx1 = dx -((5.667671E-05*dx**5+9.521310E-04*dx**4-4.787537E-01*dx**3-4.868640E+00*dx**2+1.852407E+03*dx)*math.pi*10**-9-dV)/((2.833836E-04*dx**4+3.808524E-03*dx**3-1.436261E+00*dx**2-9.737279E+00*dx+1.852407E+03)*math.pi*10**-9)
            if abs(dx1-dx) < 0.0000001:
                dx=dx1
                break
            dx=dx1
        dx = dx/1000
        x=0.24+dx
    elif V>0:
        #print "出口部"
        dV=V
        dx=dV/(0.01075**2*math.pi)
        x=0.306-dx
    else:
        print "水体積がマイナス"
        return -1
    return x

#水面位置から水面の面積を求める
def calc_S(x):
    S=0.0
    r=0.0
    dx=0.0
    if x>0.306:
        print "水面位置がボトルの外"
        return -1
    elif x<0.0:
        print "水面位置がマイナス"
        return -1
    elif x<0.065:
        print "水容量オーバー"
        return -1
    elif x<0.187:
        #print "円筒部"
        r=0.045
    elif x<0.200:
        #print "拡大部"
        dx=x-0.187
        r=0.045+dx*2.0/13.0
    elif x<0.240:
        #print "縮小部"
        dx=x-0.200
        r=0.047-dx*1.0/10.0
    elif x<0.280:
        #print "2次関数部"
        dx=(x-0.240)
        r=(-0.016834*(x*1000.0)**2+7.9672*(x*1000.0)-899.45)/1000
    elif x>=0.280:
        #print "出口部"
        r=0.01075
    S=r**2*math.pi
    return S


#----- 初期化 -----
x0 = calc_x(V0) #初期水面位置[m]
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
def fV(V):
    x=calc_x(V)
    S=calc_S(x)
    dVdt = St* math.sqrt(2/rho * (P0*(Vair0/(Vpet-V))**gamma -PA))
    return dVdt

# ルンゲクッタでVを求める
def calc_V(V):
    k=[0.0,0.0,0.0,0.0]
    k[0] = fV(V)
    k[1] = fV(V+dt*0.5*k[0])
    k[2] = fV(V+dt*0.5*k[1])
    k[3] = fV(V+dt*k[2])
    V1 = V-dt/6.0 *(k[0]+2.0*k[1]+2.0*k[2]+k[3])
    return V1


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
print "時間(t) 水面位置(x) 圧力(P) 噴出速度(v) 水面面積(S) ロケット質量(M) 推力(F) 水体積(V)"
#初期(t=0)
v=math.sqrt(2*(P-PA)/rho)
F=rho*St*v**2
print str(t),
print str(x),
print str(P),
print str(v),
print str(S),
print str(M),
print str(F),
print str(V)
t+=dt
while True:
    V1=calc_V(V)
    if V1<0:
        # ここで残りの水噴射の計算を実施
        if V>Vp[0]:
            # 水面が出口部に達していなければエラーで終了
            print "時間当たりの水体積変化が大きすぎます。時間の刻みを小さくしてください"
            break
        else:
            # 残りは噴出速度vを一定として排出時間とタンク内圧力変化を計算。
            x=Xpet[5]/1000
            Vair=Vpet
            P=P0*(Vair0/Vair)**gamma
            if P-PA<0.1:
                print "タンク内圧が大気圧まで下がりました。初期圧力を上げてください"
                break
            S=calc_S(x)
            M=m
            v=math.sqrt(2*(P-PA)/rho)
            h=V/(v*St)
            t+=h
            F=rho*St*v**2
            V=0
            print str(t),
            print str(x),
            print str(P),
            print str(v),
            print str(S),
            print str(M),
            print str(F),
            print str(V)
            print "水噴出ステージ正常終了"
            break
    V=V1
    x=calc_x(V)
    Vair=Vpet-V
    P=P0*(Vair0/Vair)**gamma
    if P-PA<0.1:
        print "タンク内圧が途中で大気圧まで下がりました。初期圧を上げてください"
        break
    v=math.sqrt(2*(P-PA)/rho)
    S=calc_S(x)
    M=m+V*rho
    F=rho*St*v**2
    print str(t),
    print str(x),
    print str(P),
    print str(v),
    print str(S),
    print str(M),
    print str(F),
    print str(V)
    t+=dt

