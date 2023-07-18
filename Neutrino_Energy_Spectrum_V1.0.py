import matplotlib.pyplot as plt
import numpy as np

# define global variables neutrino energy E/Mev, baseline distance l/km, 4 distributions
Energy = []
F_L_E_1, F_L_E_2, F_L_E_3, F_L_E_4 = [],[],[],[]
L_E = np.arange(0, 50, 0.3)
l = 58

# define the scattering cross section km^2
def sigma_E(k):
    E_e = k - (939.565-938)
    P_e = np.sqrt(E_e**2 - 0.511**2)
    sigma_E = 0.0952*(E_e*P_e)*(10**(-42))
    return sigma_E

# define reactor flux
def phi(E):
    phi_E = 0.58*np.exp(0.870 - 0.160*E - 0.091*(E**2)) + 0.30*np.exp(0.896 - 0.239*E - 0.0981*(E**2)) + 0.07*np.exp(0.976 - 0.162*E - 0.0790*(E**2)) + 0.05*np.exp(0.793 - 0.080*E - 0.1085*(E**2))
    return phi_E

# define ev^2
delta_m_21 = 7.41*(10**(-2))
delta_m_31_NO = 2.507*(10**(0))
delta_m_32_IO = -2.486*(10**(0))

# define sin^2(\theta_12)
sin_theta_12_2_NO = 0.303
sin_theta_12_2_IO = 0.303

# define sin^2(\theta_13)
sin_theta_13_2_NO = 0.02225
sin_theta_13_2_IO = 0.02223

# define sin^2(\theta_23)
sin_theta_23_2_NO = 0.451
sin_theta_23_2_IO = 0.569

P_ee_none = 1

def P_ee_only21(L_E):
    P_ee_only21 = 1-((1-sin_theta_13_2_NO)**2)*(4*sin_theta_12_2_NO*(1-sin_theta_12_2_NO))*((np.sin(delta_m_21*L_E*1.2))**2)
    return P_ee_only21

def P_ee_NO(L_E):
    P_ee_NO = 1-((1-sin_theta_13_2_NO)**2)*(4*sin_theta_12_2_NO*(1-sin_theta_12_2_NO))*((np.sin(delta_m_21*L_E*1.2))**2)-(1-sin_theta_12_2_NO)*(4*sin_theta_13_2_NO*(1-sin_theta_13_2_NO))*((np.sin(delta_m_31_NO*L_E*1.2))**2)-(sin_theta_12_2_NO)*(4*sin_theta_13_2_NO*(1-sin_theta_13_2_NO))*((np.sin((delta_m_31_NO-delta_m_21)*L_E*1.2))**2)
    return P_ee_NO

def P_ee_IO(L_E):
    P_ee_IO = 1-((1-sin_theta_13_2_IO)**2)*(4*sin_theta_12_2_IO*(1-sin_theta_12_2_IO))*((np.sin(delta_m_21*L_E*1.2))**2)-(1-sin_theta_12_2_IO)*(4*sin_theta_13_2_IO*(1-sin_theta_13_2_IO))*((np.sin((delta_m_32_IO+delta_m_21)*L_E*1.2))**2)-(sin_theta_12_2_IO)*(4*sin_theta_13_2_IO*(1-sin_theta_13_2_IO))*((np.sin(delta_m_32_IO*L_E*1.2))**2)
    return P_ee_IO

# get e
for l_e in L_E:
    e = l/l_e
    Energy.append(e)
    print(P_ee_only21(l_e))

# find 4 distributions
for l_ee, e in zip(L_E, Energy):
    # write the chi-square value
    F_L_E_1_vv = P_ee_none*phi(e)*sigma_E(e)
    F_L_E_1.append(F_L_E_1_vv)

    F_L_E_2_vv = P_ee_only21(l_ee)*phi(e)*sigma_E(e)
    F_L_E_2.append(F_L_E_2_vv)

    F_L_E_3_vv = P_ee_NO(l_ee)*phi(e)*sigma_E(e)
    F_L_E_3.append(F_L_E_3_vv)

    F_L_E_4_vv = P_ee_IO(l_ee)*phi(e)*sigma_E(e)
    F_L_E_4.append(F_L_E_4_vv)

fig, ax = plt.subplots()
ax.set_xlim([0, 10])
ax.set_ylim([0, 2*10**(-43)])
ax.plot(Energy, F_L_E_1, linestyle = '-', color = 'black', label = 'NO Oscillation')
ax.plot(Energy, F_L_E_2, linestyle = '--', color = 'b', label = '1-'+r'$\rm{P_{21}}$'+'Oscillation')
ax.plot(Energy, F_L_E_3, linestyle = '-.', color = 'g', label = 'NO')
ax.plot(Energy, F_L_E_4, linestyle = ':', color = 'r', label = 'IO')
ax.set_xlabel('L=58Km/E(Km/MeV)')
ax.set_ylabel('F(L/E)')
ax.set_title('neutrino spectrum')
plt.legend()
plt.show()
