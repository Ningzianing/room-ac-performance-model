from scipy.optimize import curve_fit
import sympy as sp
import CoolProp.CoolProp as HAP
from scipy.optimize import fsolve
import pandas as pd
def T_evp(Troom,L,BP,Q_in,):
    M_in = Q_in / 3600 * 353.25/(Troom+273.15)
    M_evp = (1 - BP) * M_in
    T_evp = Troom + L/(M_evp * 1000)
    return T_evp+273.15

def T_cnd(Tout,RHout,L,BP,Q_out,P):
    M_out = Q_out / 3600 * (353.25/(Tout+273.15))
    M_cnd = (1 - BP) * M_out
    T_cnd = Tout - (L - P)/(M_cnd * 1000)
    Q_EXT_latent = 0
    dew_point_temperature  = HAP.HAPropsSI('Tdp','T', Tout+273.15,'P', 101325,'R', RHout / 100.0)-273.15

    if (T_cnd < dew_point_temperature) & (T_cnd>0):
        h_cnd = HAP.HAPropsSI('H', 'T', Tout + 273.15, 'P', 101325, 'R', RHout / 100.0) - (L - P) / M_cnd
        T_cnd = HAP.HAPropsSI('T', 'H', h_cnd, 'P', 101325, 'R', 1) - 273.15

    if (T_cnd < dew_point_temperature) & (T_cnd<0):
        X_out = HAP.HAPropsSI('W', 'T', Tout+273.15,'P', 101325,'R', RHout / 100.0)  #
        # 求解 delta_t
        delta_t_solution = fsolve(
            lambda delta_t: X_out - (1.006 * delta_t / (1.805 + 2501 + 333.55))
                            - HAP.HAPropsSI('W', 'T', T_cnd + delta_t + 273.15, 'P', 101325, 'R', 1),
            0  # 初始猜测值
        )[0]

        T_cnd = T_cnd + delta_t_solution
        Q_EXT_latent = (M_cnd * 1000) * delta_t_solution

    return T_cnd+273.15,Q_EXT_latent

def quadratic_function(x, a, b, c):
    return a * x ** 2 + b * x + c

def calculate_cop_r(R,L,T_e,T_c,P_c,):
    cop_r = R * L * (T_e / (T_e - T_c)) / (L + P_c * R * (T_e / (T_e - T_c)))
    return cop_r

def find_equal_cop(cop_initial,Tout,RHout,BP,Q_out, R,L,T_e,P_c, tolerance=0.01, max_iterations=200):
    cop = cop_initial
    iteration = 0
    while iteration < max_iterations:
        if cop < 1:
            cop = 1
        P = L / cop
        T_c = T_cnd(Tout,RHout, L, BP, Q_out, P)[0]
        cop_r = calculate_cop_r(R,L,T_e,T_c,P_c,)
        if abs(cop_r - cop) < tolerance:
            break
        # 二分法
        cop = float((cop + cop_r) / 2.0)
        iteration += 1

    return cop


class heating_aircon():
    def __init__(self, BP_evp,BP_cnd,L_low,L_medium,L_high,P_low,P_medium,P_high,Q_in,Q_out,L_COLD,P_COLD): #室内機bypass,室外機bypass,最小能力,定格能力,最大能力,最小ec,定格ec,最大ec,室内機風量,室外機風量,低温能力,低温ec

        Tout,RHout,Troom,RHroom,Tout_cold,RHout_cold = 7,87,20,40,2,84
        P_c = sp.symbols('P_c')
        R_low    = ((L_low/P_low)      *L_low   /(L_low-(L_low/P_low)* P_c))    / (T_evp(Troom,L_low,BP_evp,Q_in)/(T_evp(Troom,L_low,BP_evp,Q_in)-T_cnd(Tout,RHout,L_low,BP_cnd,Q_out,P_low)[0]))
        R_medium = ((L_medium/P_medium)*L_medium/(L_medium-(L_medium/P_medium)* P_c)) / (T_evp(Troom,L_medium,BP_evp,Q_in)/(T_evp(Troom,L_medium,BP_evp,Q_in)-T_cnd(Tout,RHout,L_medium,BP_cnd,Q_out,P_medium)[0]))
        solution = sp.solve(R_low - R_medium, P_c)
        P_c = solution[-1]
        self.P_c =   P_c
        self.P_medium = P_medium
        R_low = ((L_low / P_low) * L_low / (L_low - (L_low / P_low) * P_c)) / (T_evp(Troom, L_low, BP_evp, Q_in) / (
                    T_evp(Troom, L_low, BP_evp, Q_in)-T_cnd(Tout,RHout, L_low, BP_cnd, Q_out, P_low)[0]))
        R_medium = ((L_medium / P_medium) * L_medium / (L_medium - (L_medium / P_medium) * P_c)) / (
                    T_evp(Troom, L_medium, BP_evp, Q_in) / (
                        T_evp(Troom, L_medium, BP_evp, Q_in) - T_cnd(Tout,RHout, L_medium, BP_cnd, Q_out, P_medium)[0]))
        R_high   = ((L_high/P_high)    *L_high  /(L_high-(L_high/P_high)* P_c))       / (T_evp(Troom,L_high,BP_evp,Q_in  )/(T_evp(Troom,L_high,BP_evp,Q_in  )-T_cnd(Tout,RHout,L_high,BP_cnd,Q_out,P_high)[0]))

        # 输入三个点的坐标
        x_data = [L_low,L_medium,L_high]
        y_data = [R_low,R_medium,R_high]
        # 进行二次拟合
        params, covariance = curve_fit(quadratic_function, x_data, y_data)

        # 得到二次拟合的系数
        self.a, self.b, self.c = params
        self.BP_evp = BP_evp
        self.Q_in = Q_in
        self.BP_cnd = BP_cnd
        self.Q_out = Q_out

        ##defrost
        converged_1 = False
        self.D = 1
        while not converged_1:
            Q_input_heat = L_COLD # inital value
            Q_average,P_ht_av = self.defrost_Q_P_cal(Tout_cold, RHout_cold, Troom, Q_input_heat, self.Q_in)

            if abs(P_ht_av - P_COLD) < 10:
                converged_1 = True
            else:
                self.D = self.D + 0.001 * (P_ht_av - P_COLD) / abs(P_ht_av - P_COLD)
        #print(self.D)

    def cal_heating_cop(self,Tout,RHout,Troom,L,Q):
        R = quadratic_function(L,self.a,self.b,self.c)
        T_e = T_evp(Troom,L,self.BP_evp,Q)
        cop_initial = 0.1
        result_cop = find_equal_cop(cop_initial,Tout,RHout,self.BP_cnd,self.Q_out, R,L,T_e,self.P_c, tolerance=0.01, max_iterations=20)
        if result_cop < 1:
            result_cop = 1

        T_c,T_c_latent = T_cnd(Tout, RHout, L, self.BP_cnd, self.Q_out, L/result_cop)
        if T_c_latent >0:
            R = quadratic_function(L,self.a,self.b,self.c) * self.D
            T_e = T_evp(Troom, L, self.BP_evp, Q)
            cop_initial = 0.1
            result_cop = find_equal_cop(cop_initial, Tout, RHout, self.BP_cnd, self.Q_out, R, L, T_e, self.P_c,
                                        tolerance=0.01, max_iterations=20)
            if result_cop < 1:
                result_cop = 1
        return result_cop

    def judge_defrost(self,Tout,RHout,Troom,L,Q):
        R = quadratic_function(L, self.a, self.b, self.c)
        T_e = T_evp(Troom, L, self.BP_evp, Q)
        cop_initial = 0.1
        result_cop = find_equal_cop(cop_initial, Tout, RHout, self.BP_cnd, self.Q_out, R, L, T_e, self.P_c,
                                    tolerance=0.01, max_iterations=20)
        if result_cop < 1:
            result_cop = 1

        T_c, T_c_latent = T_cnd(Tout, RHout, L, self.BP_cnd, self.Q_out, L / result_cop)
        if T_c_latent > 0:
            defrost_run = 1
        else:
            defrost_run = 0
        return defrost_run

    def defrost_Q_P_cal(self,Tout_cold, RHout_cold, Troom, Q_input_heat, Q_in,):
        L_COLD = Q_input_heat
        converged = False
        while not converged:
            cop_ = self.cal_heating_cop(Tout_cold, RHout_cold, Troom, Q_input_heat, Q_in)
            P_ht = Q_input_heat / cop_
            t_cnd, Q_EXT_latent = T_cnd(Tout_cold, RHout_cold, Q_input_heat, self.BP_cnd, self.Q_out, P_ht)
            r_DFtoHt = Q_EXT_latent * 333.55 / (self.P_medium * (2501 + 333.55))
            Q_input_average = Q_input_heat / (1 + r_DFtoHt)
            if Q_input_average - L_COLD < 0:
                Q_input_heat = Q_input_heat + 100
            else:
                converged = True

        P_ht_av = (P_ht  + r_DFtoHt * self.P_medium) / (1 + r_DFtoHt)
        return Q_input_average,P_ht_av

    def cal_all_conditions_cop(self,Tout,RHout,Troom,L,Q):
        defrost_run = self.judge_defrost(Tout,RHout,Troom,L,Q)
        if defrost_run == 0:
            cop = self.cal_heating_cop(Tout,RHout,Troom,L,Q)
        if defrost_run == 1:
            Q_input_average, P_ht_av = self.defrost_Q_P_cal(Tout, RHout, Troom, L, Q,)
            cop = Q_input_average/P_ht_av
        return cop


if __name__ == "__main__":
    ac = heating_aircon(0.15,0.25,700,7100,12700,100,1530,4000,1490,2530,9300,3590)
    cop = ac.cal_heating_cop(7,87,20,12700,1490) #室外温度、室外湿度、室内温度、実際負荷、実際風量
    print(cop)














