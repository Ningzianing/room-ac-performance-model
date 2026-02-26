import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import CoolProp.CoolProp as HAP
import sympy as sp

def T_evp(Troom,RHroom,L_total,BP,Q_in,): #蒸発温度計算
    dew_point_temperature  = HAP.HAPropsSI('Tdp','T', Troom+273.15,'P', 101325,'R', RHroom / 100.0)-273.15
    M_in = Q_in / 3600 * 353.25/(Troom+273.15)
    M_evp = (1 - BP) * M_in
    T_evp = Troom - L_total/(M_evp * 1005)
    if T_evp < dew_point_temperature:
        h_evp = HAP.HAPropsSI('H','T', Troom+273.15,'P', 101325,'R', RHroom / 100.0) - L_total / M_evp
        T_evp = HAP.HAPropsSI('T','H', h_evp,'P', 101325,'R', 1) -273.15
    return T_evp+273.15
def T_cnd(Tout,L,BP,Q_out,P):# 凝縮温度計算
    M_out = Q_out / 3600 * 353.25/(Tout+273.15)
    M_cnd = (1 - BP) * M_out
    T_cnd = Tout + (L + P)/(M_cnd * 1005)
    return T_cnd+273.15
def quadratic_function(x, a, b, c):
    return a * x ** 2 + b * x + c

def calculate_cop_r(R,L,T_e,T_c,P_c,):
    cop_r = R * L * (T_e / (T_c - T_e)) / (L + P_c * R * (T_e / (T_c - T_e)))
    return cop_r

def find_equal_cop(cop_initial,Tout,BP,Q_out, R,L_total,T_e,P_c, tolerance=0.01, max_iterations=200):
    cop = cop_initial
    iteration = 0
    while iteration < max_iterations:
        P = L_total / cop
        T_c = T_cnd(Tout, L_total, BP, Q_out, P)
        cop_r = calculate_cop_r(R,L_total,T_e,T_c,P_c,)
        if abs(cop_r - cop) < tolerance:
            break
        cop = (cop + cop_r) / 2.0
        iteration += 1
        #print(iteration)

    return cop


class cooling_aircon():
    def __init__(self, BP_evp,BP_cnd,L_low,L_medium,L_high,P_low,P_medium,P_high,Q_in,Q_out,):#室内機bypass　factor , 室外機bypass factor , 最小能力、定格能力、最大能力,最小EC,定格EC,最大EC,室内機風量、室外機風量
        Tout,Troom,RHroom = 35,27,47#JIS条件の室内外条件
        P_c = sp.symbols('P_c')
        R_low = ((L_low / P_low) * L_low / (L_low - (L_low / P_low) * P_c)) / (
                    T_evp(Troom, RHroom, L_low, BP_evp, Q_in) / (
                        T_cnd(Tout, L_low, BP_cnd, Q_out, P_low) - T_evp(Troom, RHroom, L_low, BP_evp, Q_in)))
        R_medium = ((L_medium / P_medium) * L_medium / (L_medium - (L_medium / P_medium) * P_c)) / (
                    T_evp(Troom, RHroom, L_medium, BP_evp, Q_in) / (
                        T_cnd(Tout, L_medium, BP_cnd, Q_out, P_medium) - T_evp(Troom, RHroom, L_medium, BP_evp, Q_in)))
        solution = sp.solve(R_low - R_medium, P_c)
        P_c = solution[-1]
        self.P_c =   P_c
        R_low = ((L_low / P_low) * L_low / (L_low - (L_low / P_low) * P_c)) / (
                    T_evp(Troom, RHroom, L_low, BP_evp, Q_in) / (
                    T_cnd(Tout, L_low, BP_cnd, Q_out, P_low) - T_evp(Troom, RHroom, L_low, BP_evp, Q_in)))
        R_medium = ((L_medium / P_medium) * L_medium / (L_medium - (L_medium / P_medium) * P_c)) / (
                T_evp(Troom, RHroom, L_medium, BP_evp, Q_in) / (
                T_cnd(Tout, L_medium, BP_cnd, Q_out, P_medium) - T_evp(Troom, RHroom, L_medium, BP_evp, Q_in)))
        R_high = ((L_high / P_high) * L_high / (L_high - (L_high / P_high) * P_c)) / (
                T_evp(Troom, RHroom, L_high, BP_evp, Q_in) / (
                T_cnd(Tout, L_high, BP_cnd, Q_out, P_high) - T_evp(Troom, RHroom, L_high, BP_evp, Q_in)))
        # 最小定格最大3点で近似する
        x_data = [L_low,L_medium,L_high]
        y_data = [R_low,R_medium,R_high]
        # 2次曲線を作る
        params, covariance = curve_fit(quadratic_function, x_data, y_data)

        # 得到二次拟合的系数
        self.a, self.b, self.c = params
        self.BP_evp = BP_evp
        self.Q_in = Q_in
        self.BP_cnd = BP_cnd
        self.Q_out = Q_out


    def cal_cooling_cop(self,Tout,Troom,RHroom,L_total,Q):#冷房負荷を計算　現在の外気温度、室内温度、室内湿度、全熱負荷により冷房COPを計算する,風量
        R = quadratic_function(L_total,self.a,self.b,self.c)
        T_e = T_evp(Troom,RHroom,L_total,self.BP_evp,Q)
        cop_initial = 0.1
        result_cop = find_equal_cop(cop_initial,Tout,self.BP_cnd,self.Q_out, R,L_total,T_e,self.P_c, tolerance=0.01, max_iterations=20)
        if result_cop < 1:
            result_cop = 1

        return result_cop

    def calculate_latent_load(self,Troom,RHroom,L_sensible,Q):
        dew_point_temperature = HAP.HAPropsSI('D', 'R', RHroom / 100, 'P', 101325, 'T',
                                              Troom + 273.15) - 273.15

        M_in = Q * 353.25/(Troom+273.15)
        M_evp = (1 - self.BP_evp) * M_in
        T_evp = Troom - 3*L_sensible/(M_evp)

        if T_evp < dew_point_temperature:
            ah_evp = HAP.HAPropsSI('W', 'T', T_evp + 273.15, 'P', 101325, 'R', 1) * 1000
        else:
            ah_evp = HAP.HAPropsSI('W', 'T', T_evp + 273.15, 'P', 101325, 'R', RHroom / 100) * 1000
        ah_room = HAP.HAPropsSI('W', 'T', Troom + 273.15, 'P', 101325, 'R', RHroom / 100) * 1000

        if T_evp > dew_point_temperature:
            latent_load = 0
        else:
            latent_load = (1 - self.BP_evp) * (ah_room - ah_evp) * Q * 1.2
        # print(t_evp,dew_point_temperature,latent_load)
        return latent_load

    def calculate_SHF(self,Troom,RHroom,L_total,Q):
        dew_point_temperature = HAP.HAPropsSI('D', 'R', RHroom / 100, 'P', 101325, 'T',
                                              Troom + 273.15) - 273.15

        T_e = T_evp(Troom,RHroom,L_total,self.BP_evp,Q,)-273.15

        if T_e < dew_point_temperature:
            ah_evp = HAP.HAPropsSI('W', 'T', T_e + 273.15, 'P', 101325, 'R', 1) * 1000
        else:
            ah_evp = HAP.HAPropsSI('W', 'T', T_e + 273.15, 'P', 101325, 'R', RHroom / 100) * 1000
        ah_room = HAP.HAPropsSI('W', 'T', Troom + 273.15, 'P', 101325, 'R', RHroom / 100) * 1000

        if T_e > dew_point_temperature:
            latent_load = 0
        else:
            latent_load = (1 - self.BP_evp) * (ah_room - ah_evp) * Q * 1.2
        SHF = (L_total-latent_load)/L_total
        return SHF




if __name__ == "__main__":
    ac = cooling_aircon(0.15,0.25,500,2200,3300,115,425,960,906,2160) #室内機bypass　factor , 室外機bypass factor , 最小能力、定格能力、最大能力,最小EC,定格EC,最大EC,室内機風量、室外機風量
    cop = ac.cal_cooling_cop(35,27,47,3400,906)                                                                 #外気温度、室内温度、室内湿度、実際負荷、実際風量

    print(cop)














