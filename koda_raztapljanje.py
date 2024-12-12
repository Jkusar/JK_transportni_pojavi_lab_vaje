import numpy as np
from tqdm import tqdm
import math
import os
import glob
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 500


plt.tight_layout()

# vhodni podatki
m0 = 0.012 # [kg] začetna masa
mD = 0.800  # [kg] masa topila
rho_s = 1751  # [kg/m^3] gostota topljenca
n_ = 1 #število delcev (privzeto 1)

# količine, ki jih moramo določiti s fitanjem
w_sat = 0.47  # [1] masni delež topljenca v topilu ob nasičenju
kw = 0.009  # [kg/m^2s] prestopnost snovi

# parametri numeričnega izračuna in arrayi, v katere zapisujemo
dt = 1  # [s] časovni korak
t_stop = 3600 * 3  # [s] čas simulacije
t_arr = [0]  # array časov v katerih izračunamo vrednosti
mt = [0 for j in range(0, int(t_stop / dt + 0.5))]  # array izračunanih vrednosti mase topljenca
i = 0  # iterator

# začetni pogoji za numerični izračun
t = 0
m = m0

# numerični izračun
for t_step in range(0, int(t_stop / dt + 0.5)):
    m = (-kw * ((36 * math.pi * (m ** 2) / (rho_s ** 2)) ** (1 / 3)) * (w_sat - (m0 - m) / mD) + m / dt) * dt
    mt[i] = m*10**3

    if m < 0.1E-3:
        break

    t += dt
    t_arr.append(t)
    i += 1

mt1 = mt[0:len(t_arr)]  # skrajšamo array mas, ker se simulacija ustavi, ko je m < 0,1 g




def eq7(t, M, kw, w_sat, mo, MD, rho_s, n_):
    return -kw*(((36*np.pi*(M**2))/(rho_s**2))**(1/3)) * (w_sat - ((n_*(mo - M))/MD))
 

def rungeKutta(M0, t0, t, dt, kw, w_sat, MD, rho_s, n_):
    n = (int)((t - t0)/dt) 
    mo = M0 #to je začetna vrednost upoštevana v enačbi in je konstantna skozi iteracije
    M = M0 #to je začetna vrednost, ki se spremeni z iteracijami
    for i in range(1, n + 1):
        k1 = dt * eq7(t0, M, kw, w_sat, mo, MD, rho_s, n_)
        k2 = dt * eq7(t0 + 0.5 * dt, M + 0.5 * k1, kw, w_sat, mo, MD, rho_s, n_)
        k3 = dt * eq7(t0 + 0.5 * dt, M + 0.5 * k2, kw, w_sat, mo, MD, rho_s, n_)
        k4 = dt * eq7(t0 + dt, M + k3, kw, w_sat, mo, MD, rho_s, n_)
 
        M = M + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        t0 = t0 + dt
    return M

#Izračun mase za določen čas
t0 = 0 #začetni čas
t = 100 # čas [s]

print ('The value of y at x is:', rungeKutta(m0, t0, t, dt, kw, w_sat, mD, rho_s, n_))
 

def get_data_fromTXT(folder_path, skip=0, convert_to_float=False):
    """
    Returns list of data from a txt file.
    
    Input parameters:
        - folder_path [string] - path to the folder with txt files. 
        If you wish to read only one txt file, specify the path to that file (it should include .txt suffix)
        - skip [int] - number of lines you wish to skip eg. txt. file has a header,... By default it is set to 0
        - convert_to_float [bool, optional] - True if you wish to convert values to floats.
    Output parameters:
        - data [dict, list] - dictionary of data or list of data in case of only one file. Keyword is file path, value is a list of data in a file
    """
    
    if ".txt" in folder_path:
        all_data =[]
        cont = open(file = folder_path, mode = "r").readlines()
        cont = cont[skip::]
        for j in cont:
            all_data.append(j.rstrip("\n"))
        return all_data      
    else:
        files = glob.glob(f"{folder_path}/**.txt") 
        all_data = {}
        for i in files:
            cont = open(file = i, mode = "r").readlines()
            cont = cont[skip::]
            data = []
            for j in cont:
                if convert_to_float:
                    number = float(j.rstrip("\n"))
                    data.append(number)
                else:
                   data.append(j.rstrip("\n"))
            #saves data from one file to a dictionary. Key is name of the .txt file
            key = os.path.basename(i).strip(".txt")
            all_data.update({key : data})
        return all_data 


#Import meritev (za primerjavo)
#Pot do txt datoteke s podatki eksperimentalnih meritev
txt_podatki = r"C:\Users\Jernej Kusar\Documents\ŠOLA 2024-2025 1. semester\Transportni pojavi\Lab. vaje\4. lab. vaja\meritve.txt"
raw_data = get_data_fromTXT(txt_podatki, skip=1)


m_exp = []
t_exp = []
for i in raw_data:
    m = i.split("\t")[2]
    if m != "/":
        t_exp.append(float(i.split("\t")[0]))
        m_exp.append(float(m))

m_rk = []
t_rk = np.arange(0, t_arr[-1], 10)
#progress_bar initialization
total_count = len(t_rk)
print("Starting Runge-Kutta calculation. This will take a while...")
progress_bar = tqdm(total = total_count, desc="Calculating")
for i in t_rk:
    # Update the progress bar
    progress_bar.update()
    
    m = rungeKutta(m0, t0, i, dt, kw, w_sat, mD, rho_s, n_)
    m_rk.append(m*10**3)

progress_bar.close()


plt.plot(np.array(t_exp)/60, m_exp, label='Meritev', marker="o")
plt.plot(np.array(t_arr)/60, mt1, label='Eulerjeva metoda')
plt.plot(np.array(t_rk)/60, m_rk, label='Runge-Kuttta metoda')
plt.xlabel('Čas [min]')
plt.ylabel('Preostala masa lizike [g]')
plt.title(f"Raztapljanje sladkorja v vodi. Časovni korak $\Delta$t = {dt} s")
plt.grid("--")
plt.legend()
plt.show()






