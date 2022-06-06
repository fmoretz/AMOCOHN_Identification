'''Function that contains the ODE system composed by equation of the AM2HN model '''


def f_Model_Deviations_AM2HN(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,y_in,k,kd,N_bac,N_S1):

    XT, X1, X2, Z, S1, S2, C = x

    if t < 20 or t > 100:
        S1in = y_in[0]      # [gCOD/L]
        S2in = y_in[1]      # [mmol/L]
        Cin  = y_in[2]      # [mmol/L]
        Zin  = y_in[3]      # [mmol/L]
        XTin = y_in[4]      # [gCOD/L]
    else:
        S1in = y_in[0]      # [gCOD/L]   
        S2in = y_in[1]      # [mmol/L]
        Cin  = y_in[2]      # [mmol/L]
        Zin  = y_in[3]      # [mmol/L]
        XTin = 1.2*y_in[4]  # [gCOD/L]

    mu1 = mu_max[0]*(S1/(S1+Ks[0]))                                                                  # Monod
    mu2 = mu_max[1]*(S2/(S2+Ks[1]+S2**2/KI2))                                                        # Haldane

    qM  = k[5]*mu2*X2                                                                                # [mmol/L/day] - Methane molar flow 
    CO2 = C + S2 - Z                                                                                 # [mmol/L]     - CO2 Dissolved
    phi = CO2 + KH*Pt + qM/kLa
    Pc  = (phi - (phi**2- 4*KH*Pt*CO2)**0.5)/(2*KH)                                                  # [atm] - Partial pressure CO2
    qC  = kLa*(CO2 - KH*Pc)                                                                          # [mmol/L/d] - Carbon Molar Flow
                                             

    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2                                                      # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                                 # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]

    return dxdt