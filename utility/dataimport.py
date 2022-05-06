# dataimport: file to import data from excel datasets
import pandas as pd

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

folder  =  r'/Users/fmoretta/Desktop/Regression/data/'
reading_path =  folder + simname + "py"+ ".xlsx"

colnames = ["HRT","S1","XT", "S2", "X1", "X2", "Z", "C","CO2","B", "pH", "q_C", "P_C", "q_CH4"]

T1 = pd.read_excel(reading_path, sheet_name = "SS_Values",header = None, names = colnames, skiprows = 1)
T2 = pd.read_excel(reading_path, sheet_name = "Influent", header = 0)
