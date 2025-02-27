import streamlit as st
import math
import pandas as pd
from datetime import datetime

st.title("Sphinx Initializer")
st.sidebar.image("Sphinxprojectlogo.png", width=150)  
menu = st.selectbox("Menu", ["Current Workspace", "View Saved Data"])
def calculate(F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox):
    Ve = (2 * T1 * (8314 / M) * (k / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    Isp = Ve / 9.81
    m_prop = F / Ve
    m_fuel = m_prop / (1 + OF)
    m_ox = OF * m_fuel
    A = math.pi * (din ** 2 / 4)
    N = m_ox / (Cd * A * (2 * (0.8 - P1) * 10**6 * D_ox) ** 0.5)
    G_ox = (4 * m_ox) / (math.pi * ID**2)  
    r = a * (G_ox ** n) 
    Cf = ((2 * k ** 2 / (k - 1)) * (2 / (k + 1)) ** ((k + 1) / (k - 1)) * (1 - (0.101325 / P1) ** ((k - 1) / k))) ** 0.5
    At = F / (Cf * P1 * 10**6)
    D_th = (4 * At / math.pi) ** 0.5
    L = m_fuel * 1000 / (D * math.pi * ID * r)
    vol_f = (m_fuel * t) / D_ox
    OD = math.sqrt((4 * vol_f / (math.pi * L)) + ID**2) * 1000  
    return {
        "Mass Flow Rate (Oxidizer)": f"{m_ox:.4f} kg/s",
        "Mass Flow Rate (Fuel)": f"{m_fuel:.4f} kg/s",
        "Number of Injectors": f"{N:.2f}",
        "Oxidizer Flux": f"{G_ox:.2f} kg/m²·s",
        "Regression Rate": f"{r:.2f} mm/s",
        "Grain Length": f"{L * 1000:.2f} mm",
        "Throat Diameter": f"{D_th * 1000:.2f} mm",
        "Exit Velocity": f"{Ve:.2f} m/s",
        "Specific Impulse": f"{Isp:.2f} s",
        "Grain Outer Diameter": f"{OD:.2f} mm"
    }

if menu == "Current Workspace":
    st.header("Current Calculation")
    st.sidebar.header("Input Parameters")
    col1, col2 = st.sidebar.columns(2)
    with col1:
        F = st.number_input("Desired Thrust in N:", value=500.0)
        P1 = st.number_input("Chamber Pressure (MPa):", value=0.4)
        T1 = st.number_input("Temperature in K:", value=3030.0)
        OF = st.number_input("O/F ratio:", value=5)
        D = st.number_input("Density (kg/m³):", value=900)
        k = st.number_input("k value (Cp/Cv):", value=1.26)
    with col2:
        a = st.number_input("Burn Rate coefficient (a):", value=0.0304)
        n = st.number_input("Pressure exponent (n):", value=0.681)
        ID = st.number_input("Grain ID (mm):", value=36) / 1000  # Convert mm to meters
        t = st.number_input("Burn time (s):", value=10.0)
        din = st.number_input("Injector Diameter (mm):", value=2.0) / 1000  # Convert mm to meters
        Cd = st.number_input("Cd:", value=0.65)
        M = st.number_input(" Molecular weight:", value=24.685)
        D_ox = st.number_input("Oxidizer density (kg/m³):", value=800.0)
    results = calculate(F, P1, T1, OF, D, k, a, n, ID, t, din, Cd, M, D_ox)

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
