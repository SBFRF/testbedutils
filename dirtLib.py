import math
import numpy as np

def dens_fun(t, s):
    # see Brad's denfun.m
    # t - temperature in degrees C
    # s - salinity in ppt
    # rho - density in kg/m3
    if type(t) == int or type(t) == float and type(s) == int or type(s) == float:
        CL = (s - 0.03) / 1.805
        if CL < 0:
            CL = 0
        else:
            pass
        rho = 1000 + 1.455 * CL - 6.5 * math.pow(10, -3) * np.power(t - 4 + 0.4 * CL, 2)
    else:
        assert t.shape == s.shape, 'Temperature and salinity arrays are not the same shape'
        CL = (s - 0.03) / 1.805
        CL[CL < 0] = 0
        rho = 1000 + 1.455 * CL - 6.5 * math.pow(10, -3) * np.power(t - 4 + 0.4 * CL, 2)
    return rho


# kinematic viscosity
def kvis_fun(t):
    # see Brad's kvisfun.m
    # t - temperature degrees C
    # kvis - kinematic viscosity in m2/s
    kvis = 1.0 * math.pow(10, -6) * (1.14 - 0.031 * (t - 15) + 6.8 * math.pow(10, -4) * np.power(t - 15, 2))
    return kvis


# fall velocity
def vfall(d50, t, s, sg):
    # see Brad's vfall.m
    # d50 is median grain size in mm
    # T is water temp in degrees C
    # S is salinity in ppt
    # sg is specific gravity of the sand grain - I added this becuase Brad hard codes it to 2.65
    g = 9.81  # acceleration due to gravity (m2/s)
    rho = dens_fun(t, s)  # get the density (kg/m3)
    kvis = kvis_fun(t)  # get the kinematic viscosity (m2/s)
    rho_s = 1000 * sg  # convert the specific gravity of sand to density (kg/m3)
    d = d50 * (1.0 / 1000)  # convert mm to m
    s = rho_s / float(rho)
    D = d * math.pow((g * (s - 1) / (math.pow(kvis, 2))), (1.0 / 3))
    w_f = (kvis / float(d)) * (np.sqrt(math.pow(10.36, 2) + 1.049 * math.pow(D, 3)) - 10.36)
    return w_f