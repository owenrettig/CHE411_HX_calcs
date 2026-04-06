# Revision 1 
# This has been compiled from two other python scripts, hxcalc_1 and hxcalc_2
# All water property data is found via IAPWS, this is a well regarded source 
# Air property data is found via CoolProp, this was done if air conditions need to be changed, see line line 114 and 115
# This file is valid for the 1-1 countercurrent heat exchanger in the RHIT UO lab

import numpy as np
from iapws import IAPWS97
import CoolProp.CoolProp as CP
import pandas as pd

# Pulling Excel data

def extract_cold_values(df):
    df = df.set_index('cold_side')
    return {key: df.loc[key, 'cold_value'] for key in df.index}

def extract_hot_values(df):
    df = df.set_index('hot_side')
    return {key: df.loc[key, 'hot_value'] for key in df.index}

def is_sheet_ready(df):
    if df.empty or df.dropna(how='all').empty:
        return False
    value_cols = [col for col in df.columns if 'value' in col.lower()]
    if not value_cols:
        return False
    for col in value_cols:
        if df[col].isna().any():
            return False
        if not pd.to_numeric(df[col], errors='coerce').notna().all():
            return False
    return True

# Main calculation

def hx_analysis(sheet_name):
    try:
        data = pd.read_excel("flow_data.xlsx", sheet_name=sheet_name)

        if not is_sheet_ready(data):
            print(f"Skipping '{sheet_name}' — sheet not filled out yet.")
            return None

        cvals = extract_cold_values(data)
        hvals = extract_hot_values(data)

        Tc_in     = cvals['T_in']
        Tc_out    = cvals['T_out']
        Pc_gauge  = cvals['P_gauge']
        cflowrate = cvals['flowrate'] * 6.309e-5  # GPM → m³/s

        Th_in     = hvals['T_in']
        Th_out    = hvals['T_out']
        Ph_gauge  = hvals['P_gauge']
        hflowrate = hvals['flowrate'] * 6.309e-5  # GPM → m³/s

        # Fluid properties 
        Tc_avg = (Tc_in + Tc_out) / 2 + 273.15      # K
        Pc     = (Pc_gauge + 14.696) * 0.00689476    # MPa (absolute)
        Th_avg = (Th_in + Th_out) / 2 + 273.15       # K
        Ph     = (Ph_gauge + 14.696) * 0.00689476    # MPa (absolute)

        cold_water = IAPWS97(T=Tc_avg, P=Pc)
        warm_water = IAPWS97(T=Th_avg, P=Ph)

        #  Heat exchanger geometry 
        St  = 0.3   * 0.0254   # m  — triangular pitch (tube spacing)
        Od  = 0.25  * 0.0254   # m  — tube outer diameter
        id_ = 0.206 * 0.0254   # m  — tube inner diameter
        sd  = 3     * 0.0254   # m  — shell inner diameter
        N_t = 60               #    — number of tubes
        l_t = 35    * 0.0254   # m  — tube length

        A_i = N_t * np.pi * id_ * l_t
        A_o = N_t * np.pi * Od  * l_t

        # EXPERIMENTAL  (hxcalc_1)

        Qc = cflowrate * cold_water.rho * cold_water.cp * (Tc_out - Tc_in)
        Qh = hflowrate * warm_water.rho * warm_water.cp * (Th_out - Th_in)
        Q_loss_exp = -(Qc + Qh)

        # CALCULATED  (hxcalc_2)

        #  Shell-side convection (ho) 
        De      = (4*(St/2)*(0.86*St) - np.pi*(Od**2)/8) / (np.pi*Od/2)   # equivalent diameter from Janna (p605)
        A_shell = np.pi*(sd/2)**2 - N_t*(Od/2)**2*np.pi                   # flow area
        Re_s    = hflowrate * De / (warm_water.nu * A_shell)
        Nu_s    = 0.36 * (Re_s**0.55) * warm_water.Prandt**(1/3)
        ho      = Nu_s * warm_water.k / De

        #  Tube-side convection (hi) — Sieder-Tate (laminar) 
        Re_t = 4 * cflowrate / (np.pi * cold_water.nu * id_ * N_t)
        Nu_t = 1.86 * (Re_t * cold_water.Prandt / (l_t / id_))**(1/3)
        hi   = Nu_t * cold_water.k / id_

        # Wall resistance 
        k_SS   = 16.2  # W / m·K  (304 stainless steel)
        R_wall = np.log(Od / id_) / (2 * np.pi * k_SS * l_t * N_t)

        # Overall heat-transfer coefficient (Uo based on outer area) 
        rhs  = 1/(hi*A_i) + R_wall + 1/(ho*A_o)
        Uo   = (rhs * A_o)**-1

        #  Log-mean temperature difference & predicted Q 
        deltaT1   = Th_in  - Tc_out
        deltaT2   = Th_out - Tc_in
        deltaTLM  = (deltaT1 - deltaT2) / np.log(deltaT1 / deltaT2)
        Q_calc    = Uo * A_o * deltaTLM

        # Free-convection loss from outer shell surface 
        T_inf  = 21       # °C  (ambient)
        rel_H  = 0.4      # —   (relative humidity)
        T_K    = T_inf + 273.15
        P_atm  = 101325   # Pa

        k_air   = CP.HAPropsSI("k",   "T", T_K, "P", P_atm, "R", rel_H)
        mu_air  = CP.HAPropsSI("mu",  "T", T_K, "P", P_atm, "R", rel_H)
        cp_air  = CP.HAPropsSI("cp",  "T", T_K, "P", P_atm, "R", rel_H)
        Vha_air = CP.HAPropsSI("Vha", "T", T_K, "P", P_atm, "R", rel_H)

        rho_air   = 1 / Vha_air
        alpha_air = k_air / (rho_air * cp_air)
        nu_air    = mu_air / rho_air
        Pr_air    = (mu_air * cp_air) / k_air

        g      = 9.81
        T_surf = data.loc[0, 'avg_surf_temp']          # °C
        T_film = (T_surf + T_inf) / 2                  # °C
        beta   = 1 / T_film                            # 1/K

        Ra_D       = (g * beta * (T_surf - T_inf) * sd**3) / (alpha_air * nu_air)
        Nu_air     = (0.6 + (0.387*Ra_D**(1/6)) / (1 + (0.559/Pr_air)**(9/16))**(8/27))**2
        h_air      = Nu_air * k_air / sd
        A_s        = sd * np.pi * l_t
        Q_loss_air = h_air * A_s * (T_surf - T_inf)   # W

        # OUTPUT

        bar = '+' + '-'*83 + '+'
        print(bar)
        print(f"|  Sheet: {sheet_name:<74}|")
        print(bar)

        print(f"| {'EXPERIMENTAL':^81} |")
        print(bar)
        print(f"|  Tube-side heat absorbed  (Qc):           {1000*Qc:>10.6g} W{' '*28}|")
        print(f"|  Shell-side heat released (Qh):           {-1000*Qh:>10.6g} W{' '*28}|")
        print(f"|  Heat loss (Q_loss = Qc + Qh):            {1000*Q_loss_exp:>10.4f} W{' '*28}|")
        print(bar)

        print(f"| {'CALCULATED':^81} |")
        print(bar)
        print(f"|  Shell-side h (ho):                       {ho:>10.6g} W/m²·K{' '*23}|")
        print(f"|  Tube-side  h (hi):                       {hi:>10.6g} W/m²·K{' '*23}|")
        print(f"|  Overall U  (Uo, outer area basis):       {Uo:>10.6g} W/m²·K{' '*23}|")
        print(f"|  Predicted heat transfer (Q_calc):        {Q_calc:>10.6g} W{' '*28}|")
        print(f"|  Free-convection shell loss (Q_loss_air): {Q_loss_air:>10.6g} W{' '*28}|")
        print(bar)
        print()

        return {
            'sheet':       sheet_name,
            'Qc_kW':       Qc,
            'Qh_kW':       -Qh,
            'Q_loss_exp_kW': Q_loss_exp,
            'ho':          ho,
            'hi':          hi,
            'Uo':          Uo,
            'Q_calc_W':    Q_calc,
            'Q_loss_air_W': Q_loss_air,
        }

    except Exception as e:
        print(f"Skipping '{sheet_name}' — {e}")
        return None

# Run all flow-rate sheets

results = []
for sheet in ['4gpm', '6gpm', '8gpm']:
    r = hx_analysis(sheet)
    if r:
        results.append(r)