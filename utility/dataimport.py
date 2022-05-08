# dataimport: file to import data from excel datasets
import pandas as pd
import numpy as np

print("Insert the name of the simulation:\n")
print("[1]: amoco_HN")
print("[2]: provaADM1")
simcode = int(input("--> "))

if simcode is 1:
	simname = "amoco_HN"
elif simcode is 2:
	simname = "provaADM1"
else:
	"please, retry, there's no data with that index"

print("Data from:",simname)

folder  =  r'/Users/fmoretta/Desktop/AMOCOHN Identification Procedure/data/'
reading_path =  folder + simname + "py"+ ".xlsx"

colnames = ["HRT","S1","XT", "S2", "X1", "X2", "Z", "C","CO2","B", "pH", "q_C", "P_C", "q_CH4"]

T1 = pd.read_excel(reading_path, sheet_name = "SS_Values",header = None, names = colnames, skiprows = 1)
T2 = pd.read_excel(reading_path, sheet_name = "Influent", header = 0)

# Get raw data
HRT   = T1["HRT"].tolist()
S1    = T1["S1"].tolist()  # [g/L] - COD
XT    = T1["XT"].tolist()
S2    = T1["S2"].tolist()  # [mmol/L] - VFA
X1   = T1["X1"].tolist()
X2   = T1["X2"].tolist()
Z     = T1["Z"].tolist()
C     = T1["C"].tolist()
CO2   = T1["CO2"].tolist()
B     = T1["B"].tolist()
pH    = T1["pH"].tolist()
qC   = T1["q_C"].tolist()  # [mmol/L/d]
PC   = T1["P_C"].tolist()  # [atm]
qM   = T1["q_CH4"].tolist()

print(f'\n{len(HRT)} data points imported')
print(f'HRT: {HRT[0]} - {HRT[-1]}')
print(f'S1: {S1}')
print(f'XT: {XT}')
print(f'S2: {S2}')
print(f'X1: {X1}')
print(f'X2: {X2}')
print(f'Z: {Z}')
print(f'C: {C}')
print(f'CO2: {CO2}')
print(f'B: {B}')
print(f'pH: {pH}')
print(f'qC: {qC}')
print(f'PC: {PC}')
print(f'qM: {qM}')
print('\n')

D = np.empty(len(HRT))
for i in range(0, len(HRT)):
	D[i] = 1/HRT[i]

S1in = float(T2["S1in"])    # [gCOD/L]
S2in = float(T2["S2in"])    # [mmol/L]
Cin  = float(T2["Cin"])    # [mmol/L]
XTin = float(T2["XTin"])    # [gCOD/L]