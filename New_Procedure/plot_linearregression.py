import numpy as np
import matplotlib.pyplot as plt
from Identification import*
from dataimport import*

marker = "s"
linewidth = 4
linestyle = "solid"

plt.subplots_adjust(
    left   = 0.1,
    bottom = 0.1,
    right  = 0.95,
    top    = 0.95,
    wspace = 0.6,
    hspace = 0.5
    )

# plt.tight_layout()
float_formatter = "{:.2f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})
fontsize = 15
numticks = 5

plt.figure(1)
x = mdl1.intercept_ + mdl1.coef_[0]*X11 + mdl1.coef_[1]*X12
plt.plot(Xr1, x, linewidth=linewidth, color="#728fa5")
plt.plot(Xr1, Yr1, marker=marker, linestyle="None", color="black")
# plt.xlabel("D*S1 [g/L/d]")
# plt.ylabel("S1 [g/L]")
plt.yticks(np.linspace(0, max(Yr1), numticks), fontsize=fontsize)
plt.xticks(np.linspace(0, np.amax(Xr1), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_1.svg', format='svg', dpi=1200)

plt.figure(2)
x = func2(Xr2,beta2[0],beta2[1],beta2[2],beta2[3])
plt.plot(x, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(Yr2, marker=marker, linestyle="None", color="black")
# plt.xlabel("D*S2 [g/L/d]", fontsize=fontsize)       # A
# plt.ylabel("S2 [g/L]", fontsize=fontsize)           # 
plt.yticks(np.linspace(0, max(Yr2), numticks), fontsize=fontsize)
plt.xticks(np.linspace(0, len(x), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_2.svg', format='svg', dpi=1200)

plt.figure(3)
y = mdl_hyd.intercept_ + mdl_hyd.coef_*X_hydr
plt.plot(X_hydr, y, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X_hydr, Y_hydr, marker=marker, linestyle="None", color="black")
# plt.xlabel("XT [g/L]", fontsize=fontsize)     # XT [g/L]
# plt.ylabel("D*(XTin-XT) [g/L/d]", fontsize=fontsize)     # Dil*(XTin-XT) [g/L/d]
plt.yticks(np.linspace(0, max(Y_hydr), numticks), fontsize=fontsize)
plt.xticks(np.linspace(0, max(X_hydr), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_3.svg', format='svg', dpi=1200)

plt.figure(4)
plt.plot(X3r, mdl3.intercept_ + mdl3.coef_*X3r, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X3r, Y3r, marker=marker, linestyle="None", color="black")
# plt.xlabel("SCO2-KH*PC [g/L/d]", fontsize=fontsize)     # C*1/(1+10**(pH-pKb))-KH*PC [g/L/d]
# plt.ylabel("qC [g/L/d]", fontsize=fontsize)     # qC [g/L/d]
plt.yticks(np.linspace(0, max(Y3r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X3r), max(X3r), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_4.svg', format='svg', dpi=1200)

plt.figure(5)
plt.plot(X4r, mdl4.intercept_ + mdl4.coef_*X4r, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X4r,Y4r, marker=marker, linestyle="None", color="black")
# plt.xlabel("G", fontsize=fontsize)     # -----
# plt.ylabel("H", fontsize=fontsize)     # -----
plt.yticks(np.linspace(0, max(Y4r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X4r), max(X4r), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_5.svg', format='svg', dpi=1200)

plt.figure(6)
plt.plot(X5r, mdl5.intercept_ + mdl5.coef_*X5r, linestyle=linestyle, linewidth=linewidth, color="#728fa5") 
plt.plot(X5r,Y5r, marker=marker, linestyle="None", color="black")
# plt.xlabel("I", fontsize=fontsize)     # -----
# plt.ylabel("L", fontsize=fontsize)     # -----
plt.yticks(np.linspace(0, max(Y5r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X5r), max(X5r), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_6.svg', format='svg', dpi=1200)

plt.figure(7)
plt.plot(X61, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X61, Y6r, marker=marker, linestyle="None", color="black")
# plt.ylabel("k6/k3*(D*(S2in-S2))+k6/k3*k2/k1*(D*(S1in-S1)+K_hyd*XT", fontsize=fontsize)     # k6/k3*(D*(S2,in -S2))+ k6/k3*k2/k1*(D*(S1,in -S1)+K_hyd*XT
# plt.xlabel("D*(S2in-S2)", fontsize=fontsize)     # D*(S2,in -S2)
plt.yticks(np.linspace(0, max(Y6r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X61), max(X61), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_7.svg', format='svg', dpi=1200)

plt.figure(8)
plt.plot(X62, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X62, Y6r, marker=marker, linestyle="None", color="black")
# plt.xlabel("D*(S1in-S1)+K_hyd*XT)", fontsize=fontsize)     # D*(S1,in -S1)+K_hyd*XT)
# plt.ylabel("P", fontsize=fontsize)     # -----
plt.yticks(np.linspace(0, max(Y6r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X62), max(X62), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_8.svg', format='svg', dpi=1200)

plt.figure(9)
plt.plot(X71, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X71,Y7r, marker=marker, linestyle="None", color="black")
# plt.ylabel("k4/k1*(D*(S1in-S1)+K_hyd*XT)+k5/k6*qM", fontsize=fontsize)     # k4/k1*(D*(S1,in -S1)+K_hyd*XT) + k5/k6*qM
# plt.xlabel("R", fontsize=fontsize)     # D*(S1,in -S1)+K_hyd*XT)
plt.yticks(np.linspace(0, max(Y7r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X71), max(X71), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_9.svg', format='svg', dpi=1200)

plt.figure(10)
plt.plot(X72, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X72,Y7r, marker=marker, linestyle="None", color="black")
# plt.xlabel("S", fontsize=fontsize)     # qM
# plt.ylabel("T", fontsize=fontsize)     # -----
plt.yticks(np.linspace(0, max(Y7r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X72), max(X72), numticks), fontsize=fontsize)
plt.grid(True, alpha=0.7, linestyle='--', color='black')
plt.savefig('New_Procedure/res_reg_fig/reg_res_10.svg', format='svg', dpi=1200)
plt.show()