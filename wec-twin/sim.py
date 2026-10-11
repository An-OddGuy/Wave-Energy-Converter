from fmpy import simulate_fmu
for b in (2.8, 5.6, 11.2):
    r = simulate_fmu(r"G:\Projects\Wave_Energy_Converter\WEC_Heave.fmu", stop_time=100,
                     start_values={"Bpto": b}, output=["E", "eta"],
                     output_interval=0.0005)
    print(b, r["E"][-1] / 100)   # mean power, W