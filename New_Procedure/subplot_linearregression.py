import numpy as np
import matplotlib.pyplot as plt
from Identification import*
from dataimport import*

plt.figure()

marker = "2"
linewidth = 1
linestyle = "solid"

plt.subplots_adjust(
    left   = 0.1,
    bottom = 0.1,
    right  = 0.95,
    top    = 0.95,
    wspace = 0.5,
    hspace = 0.5
    )

plt.tight_layout()
float_formatter = "{:.2f}".format
np.set_printoptions(formatter={'float_kind':float_formatter})
fontsize = 6.6
numticks = 4

# plt.subplot(5,2,1)
# plt.plot(Xr1, mdl1.intercept_ + mdl1.coef_[0]*X11 + mdl1.coef_[1]*X12, color="#728fa5")
# plt.plot(Xr1, Yr1, marker=marker, linestyle="None", color="black")
# plt.xlabel("D*S1 [g/L/d]")
# plt.ylabel("S1 [g/L]")

plt.subplot(3,3,1, title="(1)")
x = func2(Xr2,beta2[0],beta2[1],beta2[2],beta2[3])
plt.plot(x, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(Yr2, marker=marker, linestyle="None", color="black")
plt.xlabel("A", fontsize=fontsize)     # D*S2 [g/L/d]
plt.ylabel("B", fontsize=fontsize)     # S2 [g/L]
plt.yticks(np.linspace(0, max(Yr2), numticks), fontsize=fontsize)
plt.xticks(np.linspace(0, len(x), numticks), fontsize=fontsize)

plt.subplot(3,3,2, title="(2)")
y = mdl_hyd.intercept_ + mdl_hyd.coef_*X_hydr
plt.plot(X_hydr, y, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X_hydr, Y_hydr, marker=marker, linestyle="None", color="black")
plt.xlabel("C", fontsize=fontsize)     # XT [g/L]
plt.ylabel("D", fontsize=fontsize)     # Dil*(XTin-XT) [g/L/d]
plt.yticks(np.linspace(0, max(Y_hydr), numticks), fontsize=fontsize)
plt.xticks(np.linspace(0, max(X_hydr), numticks), fontsize=fontsize)

plt.subplot(3,3,3, title="(3)")
plt.plot(X3r, mdl3.intercept_ + mdl3.coef_*X3r, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X3r, Y3r, marker=marker, linestyle="None", color="black")
plt.xlabel("E", fontsize=fontsize)     # C*1/(1+10**(pH-pKb))-KH*PC [g/L/d]
plt.ylabel("F", fontsize=fontsize)     # qC [mmol/L/d]
plt.yticks(np.linspace(0, max(Y3r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X3r), max(X3r), numticks), fontsize=fontsize)

plt.subplot(3,3,4,title="(4)")
plt.plot(X4r, mdl4.intercept_ + mdl4.coef_*X4r, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X4r,Y4r, marker=marker, linestyle="None", color="black")
plt.xlabel("G", fontsize=fontsize)     # mu1*X1
plt.ylabel("H", fontsize=fontsize)     # Dil*(S1_in - S1) + k_hyd*XT
plt.yticks(np.linspace(0, max(Y4r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X4r), max(X4r), numticks), fontsize=fontsize)

plt.subplot(3,3,5,title="(5)")
plt.plot(X5r, mdl5.intercept_ + mdl5.coef_*X5r, linestyle=linestyle, linewidth=linewidth, color="#728fa5") 
plt.plot(X5r,Y5r, marker=marker, linestyle="None", color="black")
plt.xlabel("I", fontsize=fontsize)     # mu_2
plt.ylabel("L", fontsize=fontsize)     # q_M/X_2
plt.yticks(np.linspace(0, max(Y5r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X5r), max(X5r), numticks), fontsize=fontsize)

plt.subplot(3,3,6,title="(6)")
plt.plot(X61, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X61, Y6r, marker=marker, linestyle="None", color="black")
plt.ylabel("M", fontsize=fontsize)     # D*(S2,in -S2)
plt.xlabel("N", fontsize=fontsize)     # qM
plt.yticks(np.linspace(0, max(Y6r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X61), max(X61), numticks), fontsize=fontsize)

plt.subplot(3,3,7, title='(7)')
plt.plot(X62, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X62, Y6r, marker=marker, linestyle="None", color="black")
plt.xlabel("O", fontsize=fontsize)     # D*(S1,in -S1)+K_hyd*XT
plt.ylabel("P", fontsize=fontsize)     # qM
plt.yticks(np.linspace(0, max(Y6r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X62), max(X62), numticks), fontsize=fontsize)

plt.subplot(3,3,8, title ='(8)')
plt.plot(X71, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X71,Y7r, marker=marker, linestyle="None", color="black")
plt.ylabel("R", fontsize=fontsize)     # k4/k1*(D*(S1,in -S1)+K_hyd*XT) + k5/k6*qM
plt.xlabel("Q", fontsize=fontsize)     # qC - D*(Cin -C)
plt.yticks(np.linspace(0, max(Y7r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X71), max(X71), numticks), fontsize=fontsize)

plt.subplot(3,3,9, title ='(9)')
plt.plot(X72, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, linestyle=linestyle, linewidth=linewidth, color="#728fa5")
plt.plot(X72,Y7r, marker=marker, linestyle="None", color="black")
plt.xlabel("S", fontsize=fontsize)     # qM
plt.ylabel("T", fontsize=fontsize)     # -----
plt.yticks(np.linspace(0, max(Y7r), numticks), fontsize=fontsize)
plt.xticks(np.linspace(min(X72), max(X72), numticks), fontsize=fontsize)

plt.show()
# plt.savefig("subplot_w_labels.png", format="png", dpi=1200)