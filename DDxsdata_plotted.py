import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv("DD xs data.csv", sep=";")
print(df.columns.tolist())
print(df.head())
print(df.describe())

# peak of the curve
peak = df.loc[df["cross_section"].idxmax()]
print(peak)

#plotting the data
plt.plot(df["Incident_energy"], df["cross_section"])
plt.title("Cross Section vs Incident Energy")
plt.xscale("log")
plt.xlabel("Incident Energy (eV)")
plt.ylabel("Cross Section (barns)")
plt.savefig("cross_section.png", dpi=150, bbox_inches="tight")
plt.show()