import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import io
from datetime import datetime

def main():
    st.set_page_config(page_title="Rocket Engine Design Calculator", layout="wide")
    st.image("sphinxfinal.png", width=150)  # Adjust width as needed
    st.title("Hybrid Rocket Engine Design Calculator")
   
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Primary Parameters")
        F = st.number_input("Desired Thrust (N)", value=1000.0, format="%.2f")
        P1 = st.number_input("Combustion Chamber Pressure (MPa)", value=4.0, format="%.2f")
        T1 = st.number_input("Combustion Chamber Temperature (K)", value=3000.0, format="%.2f")
        OF = st.number_input("O/F Ratio", value=7.0, format="%.2f")
        D = st.number_input("Propellant Density (kg/m³)", value=920.0, format="%.2f")
        
    with col2:
        st.subheader("Engine Parameters")
        k = st.number_input("k Value (Cp/Cv)", value=1.2, format="%.2f")
        a = st.number_input("Fuel Regression Coefficient (a)", value=0.155, format="%.3f")
        n = st.number_input("Mass Flux Exponent (n)", value=0.5, format="%.2f")
        ID = st.number_input("Initial ID of the Grain (mm)", value=30.0, format="%.2f") / 1000  # Convert mm to meters
        t = st.number_input("Burn Time (seconds)", value=10.0, format="%.2f")
        
    col3, col4 = st.columns(2)
    
    with col3:
        st.subheader("Injector Parameters")
        din = st.number_input("Injector Diameter (mm)", value=1.0, format="%.2f") / 1000  # Convert mm to meters
        Cd = st.number_input("Discharge Coefficient of the Injector", value=0.7, format="%.2f")
        P_in = st.number_input("Inlet Pressure (MPa)", value=5.0, format="%.2f")
    
    with col4:
        st.subheader("Gas Properties")
        M = st.number_input("Molecular Weight of Exhaust Gases (g/mol)", value=22.0, format="%.2f")
        D_ox = st.number_input("Oxidizer Density (kg/m³)", value=1230.0, format="%.2f")
        T0 = st.number_input("Oxidizer Temperature (K)", value=224.0, format="%.2f")
        k_ox = st.number_input("Specific Heat Ratio for Oxidizer", value=1.27, format="%.2f")
        
    residual_percentage = st.slider("Residual Fuel Thickness (%)", min_value=1, max_value=30, value=10)
    
    if st.button("Calculate Engine Parameters"):
        with st.spinner("Calculating..."):
            results = calculate_engine_parameters(F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox, P_in, T0, k_ox, residual_percentage)
            
            # results
            st.header("Results")
            
            col_res1, col_res2 = st.columns(2)
            
            with col_res1:
                st.subheader("Performance Parameters")
                st.info(f"Mass Flow Rate of Oxidizer: {results['m_ox']:.4f} kg/s")
                st.info(f"Mass Flow Rate of Fuel: {results['m_fuel']:.4f} kg/s")
                st.info(f"Number of Injectors: {results['N_injectors']:.2f}")
                st.info(f"Regression Rate: {results['r']:.2f} mm/s")
                st.info(f"Exit Velocity: {results['Ve']:.2f} m/s")
                st.info(f"Specific Impulse: {results['Isp']:.2f} s")
                
            with col_res2:
                st.subheader("Dimensional Parameters")
                st.info(f"Grain Length: {results['L'] * 1000:.2f} mm")
                st.info(f"Throat Diameter: {results['D_th'] * 1000:.2f} mm")
                st.info(f"Final Port Diameter: {results['OD'] * 1000:.2f} mm")
                st.info(f"Outer Diameter with {residual_percentage}% residual fuel: {results['final_OD'] * 1000:.2f} mm")
                st.info(f"Total Fuel Mass: {results['total_fuel_mass']:.2f} kg")
            
            # plots
            st.subheader("Graphs")
            
            tab1, tab2, tab3, tab4 = st.tabs(["Regression Rate vs Time", "Regression Rate vs Mass Flux", 
                                              "Oxidizer Mass Flux vs Time", "Outer Diameter vs Time"])
            
            with tab1:
                fig1 = plt.figure(figsize=(10, 6))
                plt.plot(results['time_array'], results['reg_rate_array'], 'b-', label='Regression Rate')
                plt.title('Regression Rate vs Time')
                plt.xlabel('Time (s)')
                plt.ylabel('Regression Rate (mm/s)')
                plt.grid(True)
                plt.legend()
                st.pyplot(fig1)
                
            with tab2:
                fig2 = plt.figure(figsize=(10, 6))
                plt.plot(results['Gox_array'], results['reg_rate_array'], 'b-')
                plt.title('Regression Rate vs Mass Flux of Oxidizer')
                plt.xlabel('Oxidizer Mass Flux (kg/m²s)')
                plt.ylabel('Regression Rate (mm/s)')
                plt.grid(True)
                st.pyplot(fig2)
                
            with tab3:
                fig3 = plt.figure(figsize=(10, 6))
                plt.plot(results['time_array'], results['Gox_array'], 'b-', linewidth=1.5)
                plt.title('Oxidizer Mass Flux vs Time')
                plt.xlabel('Time (s)')
                plt.ylabel('Oxidizer Flux (kg/m²s)')
                plt.grid(True)
                st.pyplot(fig3)
                
            with tab4:
                fig4 = plt.figure(figsize=(10, 6))
                plt.plot(results['time_array'], results['OD_array'] * 1000, 'r-', label="Outer Diameter")
                plt.title("Outer Diameter vs Time")
                plt.xlabel("Time (s)")
                plt.ylabel("Outer Diameter (mm)")
                plt.grid(True)
                plt.legend()
                st.pyplot(fig4)
            
            #data exp
            st.subheader("Export Data")
            df = pd.DataFrame({
                "Time (s)": results['time_array'], 
                "Regression Rate (mm/s)": results['reg_rate_array'], 
                "Outer Diameter (mm)": results['OD_array'] * 1000,
                "Oxidizer Mass Flux (kg/m²s)": results['Gox_array']
            })
            
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download Data as CSV",
                data=csv,
                file_name=f"rocket_engine_data_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv"
            )
            
            # data head
            st.subheader("Data Sample")
            st.dataframe(df.iloc[::len(df)//10][:10])  # Show ~10 evenly spaced rows

def calculate_engine_parameters(F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox, P_in, T0, k_ox, residual_percentage):
    # Constants & Initial Calculations
    Ve = (2 * T1 * (8314 / M) * (k / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    Isp = Ve / 9.81
    m_prop = F / Ve
    m_fuel = m_prop / (1 + OF)
    m_ox = OF * m_fuel
    A = math.pi * (din ** 2 / 4)

    # Inputs for injector calculation
    P0 = P_in
    Pe = P1
    R_universal = 8.314      # J/mol·K
    M_N2O = 0.044013         # kg/mol (molar mass of N2O)
    R = R_universal / M_N2O  # Specific gas constant for N2O in J/kg·K

    # Effective area of one injector
    A_t = math.pi * (din ** 2) / 4
    A_eff = Cd * A_t

    # Pressure ratio
    Pr = Pe / P0

    # Mass flow rate per injector
    term1 = (P0 * 10**6 * A_eff) / (T0 ** 0.5)
    term2 = (2 * k_ox) / (R * (k_ox - 1))
    term3 = Pr ** (2 / k_ox)
    term4 = Pr ** ((k_ox + 1) / k_ox)
    m_dot_injector = term1 * math.sqrt(term2 * (term3 - term4))

    # Number of injectors needed
    N_injectors = m_ox / m_dot_injector

    G_ox = (4 * m_ox) / (math.pi * ID**2)
    r = a * (G_ox ** n)
    Cf = ((2 * k ** 2 / (k - 1)) * (2 / (k + 1)) ** ((k + 1) / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    At = F / (Cf * P1 * 10**6)
    D_th = (4 * At / math.pi) ** 0.5
    L = m_fuel * 1000 / (D * math.pi * ID * r)

    timesteps = 1000  # Reduced for better performance in web app
    time_step = t / timesteps

    # Arrays to store values
    reg_rate_array = np.zeros(timesteps)
    Gox_array = np.zeros(timesteps)
    ID_array = np.zeros(timesteps)
    OD_array = np.zeros(timesteps)
    time_array = np.linspace(time_step, t, timesteps)

    # Initial Outer Diameter
    r_t = ID / 2
    OD = 2 * r_t

    # Simulation loop
    for i in range(timesteps):
        G_ox = (4 * m_ox) / (math.pi * ID**2)
        reg_rate = a * G_ox ** n
        reg_rate_array[i] = reg_rate
        Gox_array[i] = G_ox
        ID += 2 * (reg_rate * 0.001) * time_step
        ID_array[i] = ID
        r_t = (a * 0.001 * (2 * n + 1) * (m_ox / math.pi) ** n * time_array[i] + (0.5 * ID) ** (2 * n + 1)) ** (1 / (2 * n + 1))
        OD = 2 * r_t
        OD_array[i] = OD

    # Calculate final outer diameter with residual
    final_port_radius = OD_array[-1] / 2
    residual_thickness = (residual_percentage / 100) * final_port_radius
    final_OD = 2 * (final_port_radius + residual_thickness)

    # Calculate the mass of fuel
    grain_radius = final_OD / 2
    grain_volume = math.pi * grain_radius**2 * L
    total_fuel_mass = D * grain_volume

    return {
        'Ve': Ve,
        'Isp': Isp,
        'm_ox': m_ox,
        'm_fuel': m_fuel,
        'N_injectors': N_injectors,
        'r': r,
        'D_th': D_th,
        'L': L,
        'OD': OD,
        'final_OD': final_OD,
        'total_fuel_mass': total_fuel_mass,
        'time_array': time_array,
        'reg_rate_array': reg_rate_array,
        'Gox_array': Gox_array,
        'OD_array': OD_array
    }

if __name__ == "__main__":
    main()
