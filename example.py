from ac_cooling import cooling_aircon
from ac_heating import heating_aircon


# Cooling model initialization
ac_cooling = cooling_aircon(
    0.15, 0.25,              # indoor unit bypass factor, outdoor unit bypass factor
    500, 2200, 3300,         # minimum capacity, rated capacity, maximum capacity [W]
    115, 425, 960,           # minimum power, rated power, maximum power [W]
    906, 2160                # indoor airflow rate, outdoor airflow rate
)

# Calculate cooling COP under actual operating condition
cop_cooling = ac_cooling.cal_cooling_cop(
    32,      # outdoor temperature [°C]
    25,      # indoor temperature [°C]
    60,      # indoor relative humidity [%]
    2000,    # actual cooling load [W]
    906      # actual airflow rate
)

# Heating model initialization
ac_heating = heating_aircon(
    0.15, 0.25,              # indoor unit bypass factor, outdoor unit bypass factor
    700, 7100, 12700,        # minimum capacity, rated capacity, maximum capacity [W]
    100, 1530, 4000,         # minimum power, rated power, maximum power [W]
    1490, 2530,              # indoor airflow rate, outdoor airflow rate
    9300, 3590               # low-temperature heating capacity, low-temperature heating power consumption
)

# Calculate heating COP under actual operating condition
cop_heating = ac_heating.cal_heating_cop(
    5,       # outdoor temperature [°C]
    70,      # outdoor relative humidity [%]
    25,      # indoor temperature [°C]
    4000,   # actual heating load [W]
    1490     # actual airflow rate
)

# Power consumption
power_cooling = 2000 / cop_cooling
power_heating = 4000 / cop_heating


print("===== Cooling Condition =====")
print("Outdoor: 32°C")
print("Indoor : 25°C, 60% RH")
print("Load   : 2000 W")
print(f"COP    : {cop_cooling:.3f}")
print(f"Power  : {power_cooling:.1f} W")

print("\n===== Heating Condition =====")
print("Outdoor: 5°C, 70% RH")
print("Indoor : 25°C")
print("Load   : 4000 W")
print(f"COP    : {cop_heating:.3f}")
print(f"Power  : {power_heating:.1f} W")


