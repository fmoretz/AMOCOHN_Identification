# Expressions to be regressed against the data
import numpy as np
from utility.PhysConstants import *
from utility.ReactorConf import *
from typing import List

def Methane(x, alpha: float, biomass_2: float, Dilution: float, Methane_flow: float) -> float:
    # Define the variables
    k6 = x
    
    # Define the expressions
    return Methane_flow / biomass_2 - alpha * k6 * Dilution

def Carbon_Dioxide(x, henry_constant: float, pH: float, Kb: float, Inorganic_Carbon: float, Carbon_Pressure: float, CO2_flow: float) -> float:
    # Define the variables
    kLa = x
    
    # Define the expressions
    return CO2_flow - kLa * ( Inorganic_Carbon / (1+10**(pH)) - henry_constant * Carbon_Pressure )

def Inorganic_Carbon(x, alpha: float, Dilution: float, biomass_1: float, biomass_2: float, IC_content: float, Initial_IC: float, CO2_flow: float) -> float:
    # Define the variables
    k4, k5 = x
    
    # Define the expressions
    return CO2_flow - Dilution * (Initial_IC - IC_content) - k4 * alpha * Dilution * biomass_1 - k5 * alpha * Dilution * biomass_2

def Hydrolysis(x, Particulate_contet: float, Initial_particulate: float, Dilution: float):
    # Define the variables
    khyd = x
    
    # Define the expressions
    return Dilution * (Initial_particulate - Particulate_contet) - khyd * Particulate_contet

def Acidogenesis(x, khyd, Substrate_1: float, Initial_S1: float, Dilution: float, alpha: float, biomass_1: float, Particulate_contet: float) -> float:
    # Define the variables
    k1 = x
    
    # Define the expressions
    return Dilution * (Initial_S1 - Substrate_1) - k1 * alpha * Dilution * biomass_1 + khyd * Particulate_contet
    
def Methanogenesis(x, Substrate_2: float, Initial_S2: float, Dilution: float, alpha: float, biomass_2: float, biomass_1: float) -> float:
    # Define the variables
    k2, k3 = x
    
    # Define the expressions
    return Dilution * (Initial_S2 - Substrate_2) + k2 * alpha * Dilution * biomass_1 - k3 * alpha * Dilution * biomass_2
    
    
    
    