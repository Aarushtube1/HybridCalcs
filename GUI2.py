import streamlit as st
import math
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

st.title("Sphinx Initializer")
st.sidebar.image("Sphinxprojectlogo.png", width=150)  
menu = st.selectbox("Menu", ["Current Workspace", "View Saved Data"])

def calculate(P_in, F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox, final_OD):
    # Constants & Initial Calculations
    Ve = (2 * T1 * (8314 / M) * (k / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    Isp = Ve / 9.81
    m_prop = F / Ve
    m_fuel = m_prop / (1 + OF)
    m_ox = OF * m_fuel 
    A = math.pi * (din ** 2 / 4)
    N = m_ox / (Cd * A * (2 * (P_in - P1) * 10**6 * D_ox) ** 0.5)
    G_ox = (4 * m_ox) / (math.pi * ID**2)
    r = a * (G_ox ** n)
    Cf = ((2 * k ** 2 / (k - 1)) * (2 / (k + 1)) ** ((k + 1) / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    At = F / (Cf * P1 * 10**6)
    D_th = (4 * At / math.pi) ** 0.5
    L = m_fuel * 1000 / (D * math.pi * ID * r)

    # Calculate Oxidizer Volume and Fuel Volume
    V_ox = m_ox / D_ox  # Oxidizer volume
    V_fuel = m_fuel / D  # Fuel volume

    # Calculate Total Propellant Mass
    m_prop_total = m_ox + m_fuel

    # Calculate Expansion Ratio 
    P_e = 0.101325  # Exit pressure in MPa 
    Ae_At = ((2 / (k + 1)) ** (1 / (k - 1))) * ((P1 / P_e) ** (1 / k)) * (((k + 1) / (k - 1)) * (1 - (P_e / P1) ** ((k - 1) / k))) ** 0.5
    Ae = Ae_At * At  # Exit area
    expansion_ratio = Ae / At  # Expansion ratio

    return {
        "Mass Flow Rate (Oxidizer)": f"{m_ox:.4f} kg/s",
        "Mass Flow Rate (Fuel)": f"{m_fuel:.4f} kg/s",
        "Number of Injectors": f"{N:.2f}",
        "Regression Rate": f"{r:.2f} mm/s",
        "Grain Length": f"{L * 1000:.2f} mm",
        "Throat Diameter": f"{D_th * 1000:.2f} mm",
        "Exit Velocity": f"{Ve:.2f} m/s",
        "Specific Impulse": f"{Isp:.2f} s",
        "Grain Outer Diameter (OD)": f"{final_OD * 1000:.2f} mm", 
        "Oxidizer Volume": f"{V_ox:.4f} m³",
        "Oxidizer Mass": f"{m_ox:.4f} kg",
        "Fuel Volume": f"{V_fuel:.4f} m³",
        "Fuel Mass": f"{m_fuel:.4f} kg",
        "Total Propellant Mass": f"{m_prop_total:.4f} kg",
        "Expansion Ratio": f"{expansion_ratio:.2f}"
    }

if menu == "Current Workspace":
    st.header("Current Calculation")
    st.sidebar.header("Input Parameters")
    col1, col2 = st.sidebar.columns(2)
    with col1:
        F = st.number_input("Desired Thrust in N:", value=250.0)
        P1 = st.number_input("Chamber Pressure (MPa):", value=0.4)
        T1 = st.number_input("Temperature in K:", value=3030.0)
        OF = st.number_input("O/F ratio:", value=5.50)
        D = st.number_input("Density (kg/m³):", value=918)
        k = st.number_input("k value (Cp/Cv):", value=1.26)
        P_in = st.number_input("Inlet Pressure (MPa):", value=0.7)
    with col2:
        a = st.number_input("Burn Rate coefficient (a):", value=0.0304)
        n = st.number_input("Pressure exponent (n):", value=0.681)
        ID = st.number_input("Grain ID (mm):", value=20) / 1000  
        t = st.number_input("Burn time (s):", value=5.0)
        din = st.number_input("Injector Diameter (mm):", value=4.0) / 1000 
        Cd = st.number_input("Cd:", value=0.62)
        M = st.number_input("Molecular weight:", value=24.880)
        D_ox = st.number_input("Oxidizer density (kg/m³):", value=14.7)

    timesteps = 1000
    burn_time = t
    time_step = burn_time / timesteps

    reg_rate_array = np.zeros(timesteps)
    Gox_array = np.zeros(timesteps)
    ID_array = np.zeros(timesteps)
    OD_array = np.zeros(timesteps)
    time_array = np.linspace(time_step, burn_time, timesteps)

    for i in range(timesteps):
        G_ox = (4 * OF * ((F / ((2 * T1 * (8314 / M) * (k / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5)) / (1 + OF))) / (math.pi * ID**2)
        reg_rate = a * G_ox ** n
        reg_rate_array[i] = reg_rate
        Gox_array[i] = G_ox
        ID += 2 * (reg_rate * 0.001) * time_step
        ID_array[i] = ID
        r_t = (a * 0.001 * (2 * n + 1) * (OF * ((F / ((2 * T1 * (8314 / M) * (k / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5)) / (1 + OF))    / math.pi) ** n * time_array[i] + (0.5 * ID) ** (2 * n + 1)) ** (1 / (2 * n + 1))
        OD = 2 * r_t
        OD_array[i] = OD

    # Final OD value from the simulation loop
    final_OD = OD_array[-1]

    # Calculate results using the final OD value
    results = calculate(P_in, F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox, final_OD)

    st.header("Output Parameters")
    output_data = {
        "Parameter": list(results.keys()),
        "Value": list(results.values())
    }
    output_df = pd.DataFrame(output_data)
    st.table(output_df)

    calculation_name = st.text_input("Name this calculation:", "Calculation")

    if st.button("Save Current Calculation"):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        saved_data = {
            "Name": calculation_name,
            "Timestamp": timestamp,
            **results
        }
        if "saved_calculations" not in st.session_state:
            st.session_state.saved_calculations = []
        st.session_state.saved_calculations.append(saved_data)
        st.success("Calculation saved!")

    st.header("Plots")

    plt.style.use('dark_background')

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 20), facecolor='#0E1117')

    line_color = '#1f77b4'
    grid_color = '#2A2A2A'
    text_color = 'white'

    # Regression Rate vs Time
    ax1.plot(time_array, reg_rate_array, color=line_color, label='Regression Rate')
    ax1.set_title('Regression Rate vs Time', color=text_color, fontsize=14, pad=10)
    ax1.set_xlabel('Time (s)', color=text_color, fontsize=12)
    ax1.set_ylabel('Regression Rate (mm/s)', color=text_color, fontsize=12)
    ax1.set_xlim([0, burn_time])
    ax1.set_ylim([0, 2])
    ax1.grid(True, color=grid_color, linestyle='--', alpha=0.7)
    ax1.legend(facecolor='#0E1117', edgecolor='white', fontsize=10)

    # Regression Rate vs Mass Flux of Oxidizer
    ax2.plot(Gox_array, reg_rate_array, color=line_color)
    ax2.set_title('Regression Rate vs Mass Flux of Oxidizer', color=text_color, fontsize=14, pad=10)
    ax2.set_xlabel('Oxidizer Mass Flux (kg/m²s)', color=text_color, fontsize=12)
    ax2.set_ylabel('Regression Rate (mm/s)', color=text_color, fontsize=12)
    ax2.grid(True, color=grid_color, linestyle='--', alpha=0.7)

    # Oxidizer Mass Flux vs Time
    ax3.plot(time_array, Gox_array, color=line_color, linewidth=1.5)
    ax3.set_title('Oxidizer Mass Flux vs Time', color=text_color, fontsize=14, pad=10)
    ax3.set_xlabel('Time (s)', color=text_color, fontsize=12)
    ax3.set_ylabel('Oxidizer Flux (kg/m²s)', color=text_color, fontsize=12)
    ax3.grid(True, color=grid_color, linestyle='--', alpha=0.7)

    # Outer Diameter vs Time
    ax4.plot(time_array, OD_array, color=line_color, label="Outer Diameter")
    ax4.set_title("Outer Diameter vs Time", color=text_color, fontsize=14, pad=10)
    ax4.set_xlabel("Time (s)", color=text_color, fontsize=12)
    ax4.set_ylabel("Outer Diameter (mm)", color=text_color, fontsize=12)
    ax4.grid(True, color=grid_color, linestyle='--', alpha=0.7)
    ax4.legend(facecolor='#0E1117', edgecolor='white', fontsize=10)

    plt.tight_layout(pad=3.0)

    st.pyplot(fig)

elif menu == "View Saved Data":
    st.header("Past Calculations")
    if "saved_calculations" in st.session_state and st.session_state.saved_calculations:
        saved_data = []
        for calculation in st.session_state.saved_calculations:
            row = {
                "Name": calculation["Name"],
                "Timestamp": calculation["Timestamp"],
                **{param: value for param, value in calculation.items() if param not in ["Name", "Timestamp"]}
            }
            saved_data.append(row)
        saved_df = pd.DataFrame(saved_data)
        st.table(saved_df)
    else:
        st.write("No past calculations available.")
