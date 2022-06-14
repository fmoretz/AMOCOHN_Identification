
from main_comparison import*
import matplotlib.pyplot as plt

gridcolor = '#c8cbdb'
gridalpha = 0.8
gridstyle = '--'
linecolor = '#2A547E' #'#'#232F3A #728fa5 366BA1 6699CC #2A547E

xaxisticks = np.linspace(0, tspan[-1], 9)

legendloc ='upper right'
legendfontsize = 6.6

fig1 = plt.figure(1)
plt.subplots_adjust(
    left   = 0.1,
    bottom = 0.1,
    right  = 0.9,
    top    = 0.9,
    wspace = 0.3,
    hspace = 0.3
    )

plt.tight_layout()

plt.subplot(2,2,1)
plt.title('AM2HN Procedure')
plt.plot(tspan,q_M_or,linestyle = 'solid', color = linecolor, label = 'qM')
plt.plot(tspan,q_C_or,linestyle = 'dashed',color = linecolor, label = 'qC')
plt.ylim(0,70)
plt.xlabel('Time [d]')
plt.ylabel('Molar Flow [mmol/L/d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks=np.linspace(0,tspan[-1],5))
plt.grid(alpha=gridalpha, color=gridcolor, linestyle=gridstyle)
plt.legend(frameon=False, fontsize='small')

plt.subplot(2,2,2)
plt.title('New Procedure')
plt.plot(tspan,q_M,linestyle = 'solid', color = linecolor, label = 'qM')
plt.plot(tspan,q_C,linestyle = 'dashed',color = linecolor, label = 'qC')
plt.xlabel('Time [d]')
plt.ylabel('Molar Flow [mmol/L/d]')
plt.ylim(0,70)
plt.xlim(0,tspan[-1])
plt.xticks(ticks=np.linspace(0,tspan[-1],5))
plt.grid(alpha=gridalpha, color=gridcolor, linestyle=gridstyle)
plt.legend(frameon=False, fontsize='small')

plt.subplot(2,2,3)
plt.plot(tspan, X1_or, linestyle = 'solid', color = linecolor, label = 'X1')
plt.plot(tspan, X2_or, linestyle = 'dashed', color = linecolor , label = 'X2')
plt.ylabel('Bacteria Concentration  [g/L]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.ylim(0,1.45)
plt.xticks(ticks=np.linspace(0,tspan[-1],5))
plt.grid(alpha=gridalpha, color=gridcolor, linestyle=gridstyle)
plt.legend( frameon=False, fontsize='small')

plt.subplot(2,2,4)
plt.plot(tspan, X1, linestyle = 'solid', color = linecolor, label = 'X1')
plt.plot(tspan, X2, linestyle = 'dashed', color = linecolor, label = 'X2')
plt.ylabel('Bacteria Concentration  [g/L]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.ylim(0,1.45)
plt.xticks(ticks=np.linspace(0,tspan[-1],5))
plt.grid(alpha=gridalpha, color=gridcolor, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize='small')

# plt.savefig("results_flows.png", format="png", dpi=1200)

fig2 = plt.figure(2)

plt.subplot(5,2,1)
plt.plot(tspan, S1_or, color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, S1, color=linecolor, label="New")
plt.ylabel('S1 [g/L]')
plt.xlabel('Time [d]')
plt.xticks(ticks = xaxisticks, labels = [])
plt.xlim(0,tspan[-1])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle, animated=True)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,2)
plt.plot(tspan, S1_or/S1_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, S1/S1[0], color=linecolor, label="New")
plt.ylabel('S1* [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,3)
plt.plot(tspan, S2_or, color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, S2, color=linecolor, label="New")
plt.ylabel('S2 [mmol/L]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels =[])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,4)
plt.plot(tspan, S2_or/S2_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, S2/S2[0], color=linecolor, label="New")
plt.ylabel('S2* [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,5)
plt.plot(tspan, CO2_or, color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, CO2, color=linecolor, label="New")
plt.ylabel('CO2 [mmol/L]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.ylim(top=CO2_or[-1]*1.5)
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,6)
plt.plot(tspan, CO2_or/CO2_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, CO2/CO2[0], color=linecolor, label="New")
plt.ylabel('CO2* [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,7)
plt.plot(tspan, C_or, color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, C, color=linecolor, label="New")
plt.ylabel('C [mmol/L]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,8)
plt.plot(tspan, C_or/C_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, C/C[0], color=linecolor, label="New")
plt.ylabel('C* [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.xticks(ticks = xaxisticks, labels = [])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,9)
plt.plot(tspan, pH_or, color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, pH, color=linecolor, label="New")
plt.ylabel('pH [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.ylim(top=8.3)
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

plt.subplot(5,2,10)
plt.plot(tspan, pH_or/pH_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, pH/pH[0], color=linecolor, label="New")
plt.ylabel('pH* [-]')
plt.xlabel('Time [d]')
plt.xlim(0,tspan[-1])
plt.grid(color=gridcolor, alpha=gridalpha, linestyle=gridstyle)
plt.legend(loc=legendloc, frameon=False, fontsize= legendfontsize)

fig2.align_ylabels()


#plt.savefig("results_SSratios.png", format="png", dpi=1200)

plt.figure(3)
plt.subplot(2,2,1)
plt.plot(tspan, q_C_or/q_C_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, q_C/q_C[0], color=linecolor, label="New")

plt.subplot(2,2,2)
plt.plot(tspan, q_M_or/q_M_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, q_M/q_M[0], color=linecolor, label="New")

plt.subplot(2,2,3)
plt.plot(tspan, X1_or/X1_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, X1/X1[0], color=linecolor, label="New")

plt.subplot(2,2,4)
plt.plot(tspan, X2_or/X2_or[0], color=linecolor, linestyle = 'dashed', label="AM2HN")
plt.plot(tspan, X2/X2[0], color=linecolor, label="New")


plt.show()