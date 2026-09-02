# propellant_trade.py
# Trading different propellants
# Author: Marlow Nedelchev

import itertools
import numpy as np
from rocketcea.cea_obj import CEA_Obj

# === FUNCTIONS ===

def print_itr(Pc: float, eps: float, MR: float):
    IspVac, Cstar, Tc, MW, gamma = C.get_IvacCstrTc_ChmMwGam(Pc=Pc, MR=MR, eps=eps)

    print( '%8.1f %8.1f   %8.1f       %8.1f      %8.1f    %8.1f %8.2f  %8.4f '%\
           (Pc, eps, MR, IspVac, Cstar, Tc, MW, gamma))

def print_table(C, Pc_arr: np.ndarray, eps_arr: np.ndarray, MR_arr: np.ndarray):
	print(' Pc(psia) AreaRatio  MixtureRatio   IspVac(sec)  Cstar(ft/sec) Tc(degR)  MolWt    gamma')
	# iteration loop and print rows
	for pc, eps, mr in itertools.product(Pc_arr, eps_arr, MR_arr):
		IspVac, Cstar, Tc, MW, gamma = C.get_IvacCstrTc_ChmMwGam(Pc=pc, MR=mr, eps=eps)
		print( '%8.1f %8.1f   %8.1f       %8.1f      %8.1f    %8.1f %8.2f  %8.4f '%\
         		(pc, eps, mr, IspVac, Cstar, Tc, MW, gamma))


# === TABLE GENERATION ===

pc_  = [100, 200, 300]
eps_ = [ 40]
mr_  = [  1,   2,   3,   4]

c_loxeth = CEA_Obj(oxName='LOX', fuelName='Ethanol')
c_loxmet = CEA_Obj(oxName='LOX', fuelName='CH4')
c_loxipa = CEA_Obj(oxName='LOX', fuelName='Isopropanol')
c_loxker = CEA_Obj(oxName='LOX', fuelName='Kerosene')

print('LOX-ETHANOL TABLE')
print_table(c_loxeth, pc_, eps_, mr_)

print('LOX-METHANE TABLE')
print_table(c_loxmet, pc_, eps_, mr_)

print('LOX-IPA TABLE')
print_table(c_loxipa, pc_, eps_, mr_)

print('LOX-KEROSENE TABLE')
print_table(c_loxker, pc_, eps_, mr_)

