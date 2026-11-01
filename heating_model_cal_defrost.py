from scipy.optimize import curve_fit
import sympy as sp
import CoolProp.CoolProp as HAP
from scipy.optimize import fsolve
import pandas as pd

COP_MAX = 20.0
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
        self.L_high = L_high
        self.P_high = P_high
        self.Vf = self.cal_Vf_heating()
        ##defrost
        converged_1 = False
        self.D = 1
        while not converged_1:
            Q_input_heat = L_COLD # inital value
            Q_input_average,P_ht_av,Q_input_heat = self.defrost_Q_P_cal(Tout_cold, RHout_cold, Troom, Q_input_heat, self.Q_in)

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
        return min(max(float(result_cop), 1.0), COP_MAX)

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
        return Q_input_average,P_ht_av,Q_input_heat

    def cal_all_conditions_cop(self,Tout,RHout,Troom,L,Q):
        defrost_run = self.judge_defrost(Tout,RHout,Troom,L,Q)
        if defrost_run == 0:
            cop = self.cal_heating_cop(Tout,RHout,Troom,L,Q)
        if defrost_run == 1:
            Q_input_average, P_ht_av,_ = self.defrost_Q_P_cal(Tout, RHout, Troom, L, Q,)
            cop = Q_input_average/P_ht_av
        return float(cop)

    def cal_Vf_heating(self):
        fluid = 'R32'
        T_inside = T_evp(20,self.L_high,self.BP_evp,self.Q_in,)
        T_outside,_ = T_cnd(7,87,self.L_high,self.BP_cnd,self.Q_out,self.P_high)

        P_inside = HAP.PropsSI('P', 'T', T_inside, 'Q', 0, fluid)
        P_outside = HAP.PropsSI('P', 'T', T_outside, 'Q', 1, fluid)
        H_compressor_in = HAP.PropsSI('H', 'T', T_outside + 5, 'P', P_outside, fluid)  # 压缩机入口
        S_compressor_in = HAP.PropsSI('S', 'T', T_outside + 5, 'P', P_outside, fluid)  # 压缩机入口熵
        H_compressor_out_s = HAP.PropsSI('H', 'P', P_inside, 'S', S_compressor_in, fluid)
        H_compressor_out_s = H_compressor_in + (H_compressor_out_s - H_compressor_in) / 0.8 # 等エントロピー効率0.8と想定
        H_expansion_in = HAP.PropsSI('H', 'T', T_inside - 5, 'P', P_inside, fluid)  # 膨胀阀入口
        density = HAP.PropsSI('D', 'H', H_compressor_in, 'P', P_outside, fluid)
        self.Vf = self.L_high / (density * (H_compressor_out_s - H_expansion_in))
        return self.Vf

    def cal_max_heating(self, Tout, Troom, RHout, Q_in, Q_out, ):
        # Return the maximum time-averaged heating capacity for the condition.
        # Frosting conditions are reduced by the defrost-to-heating time ratio.
        L_initial = self.L_high
        fluid = 'R32'
        tolerance = 10  # 容许误差
        max_iterations = 50  # 最大迭代次数
        iteration = 0

        while iteration < max_iterations:
            cop = self.cal_heating_cop(Tout,RHout,Troom,L_initial,Q_in)
            T_inside = T_evp(Troom, L_initial, self.BP_evp, Q_in, )
            T_outside, _ = T_cnd(Tout, RHout, L_initial, self.BP_cnd, Q_out, L_initial/cop)

            P_inside = HAP.PropsSI('P', 'T', T_inside, 'Q', 0, fluid)
            P_outside = HAP.PropsSI('P', 'T', T_outside, 'Q', 1, fluid)
            H_compressor_in = HAP.PropsSI('H', 'T', T_outside + 5, 'P', P_outside, fluid)  # 压缩机入口
            S_compressor_in = HAP.PropsSI('S', 'T', T_outside + 5, 'P', P_outside, fluid)  # 压缩机入口熵
            H_compressor_out_s = HAP.PropsSI('H', 'P', P_inside, 'S', S_compressor_in, fluid)
            H_compressor_out_s = H_compressor_in + (H_compressor_out_s - H_compressor_in) / 0.8

            H_expansion_in = HAP.PropsSI('H', 'T', T_inside - 5, 'P', P_inside, fluid)  # 膨胀阀入口
            density = HAP.PropsSI('D', 'H', H_compressor_in, 'P', P_outside, fluid)
            max_heating = density * self.Vf * (H_compressor_out_s - H_expansion_in)

            if abs(max_heating - L_initial) < tolerance:
                # print(f"Convergence achieved: L_initial = {L_initial}, max_cooling = {max_cooling}")
                break

            L_initial = (L_initial + max_heating) / 2  # 更新L_initial的值
            iteration += 1
        # L_initial is the maximum capacity while the unit is actively heating.
        # In a frosting condition, convert it to the average capacity over the
        # complete heating/defrost cycle using the same time-ratio model as
        # defrost_Q_P_cal().
        if self.judge_defrost(Tout, RHout, Troom, L_initial, Q_in) == 0:
            return L_initial

        cop = self.cal_heating_cop(Tout, RHout, Troom, L_initial, Q_in)
        P_ht = L_initial / cop
        _, Q_EXT_latent = T_cnd(Tout, RHout, L_initial, self.BP_cnd, Q_out, P_ht)
        r_DFtoHt = Q_EXT_latent * 333.55 / (self.P_medium * (2501 + 333.55))
        return L_initial / (1 + r_DFtoHt)

if __name__ == "__main__":
    ac = heating_aircon(0.15,0.25,700,7100,12700,100,1530,4000,1490,2530,9300,3590)
    #max_heating =ac.cal_max_heating(7, 20, 87, 1490, 2530, )
    #defrost = ac.judge_defrost(7, 87, 20,12700, 1490, )
    cop = ac.cal_heating_cop(7,87,20,12700,1490) #室外温度、室外湿度、室内温度、実際負荷、実際風量
    print(cop)
















