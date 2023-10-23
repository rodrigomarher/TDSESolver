import numpy as np
import matplotlib.pyplot as plt

w0 = 0.057/2.0

def LoadComplexData(file,**genfromtext_args):
    """
    Load complex data in the C++ format in numpy.
    """
    array_as_strings = np.genfromtxt(file,dtype=str,**genfromtext_args)
    complex_parser = np.vectorize(lambda x: complex(*eval(x)))
    return complex_parser(array_as_strings)



path = "results/"
acc_i = LoadComplexData(f"{path}/acc_i.dat")
acc_j = LoadComplexData(f"{path}/acc_j.dat")
acc_k = LoadComplexData(f"{path}/acc_k.dat")
t = np.genfromtxt(f"{path}/time.dat")
efield_k = np.genfromtxt(f"{path}/Efield_k.dat")
dt = t[1]-t[0]
nt = t.shape[0]
period = 2*np.pi/w0

esp_i = np.abs(np.fft.rfft(acc_i.real))**2
esp_j = np.abs(np.fft.rfft(acc_j.real))**2
esp_k = np.abs(np.fft.rfft(acc_k.real))**2

dw = 2*np.pi/(nt*dt)
wmax = nt*dw/2.0
w = np.arange(0,wmax+dw,dw)/w0
if(w.shape[0] != esp_k.shape[0]):
    if(w.shape[0]-1 == esp_k.shape[0]):
        w = w[:-1]
    elif w.shape[0]+1==esp_k.shape[0]:
        esp_k = esp_k[:-1]
    else:
        print("Error. Frequency shape do not match")

fig = plt.figure(figsize=(12,5))
gs = plt.GridSpec(nrows=1, ncols=2, width_ratios=[1,1])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])

ax0.plot(t/period, efield_k,linewidth=2)
ax0.set_xlabel("Time (optical cycles)", fontsize=16)
ax0.set_ylabel("Efield z (a.u.)", fontsize=16)
ax0.tick_params(which='major', axis='both', labelsize=16)
ax0.set_title("Driving electric field", fontsize=20)

ax1.plot(w, esp_k, linewidth=2)
ax1.set_yscale('log')
ax1.set_xlim(0,200)
ax1.set_ylabel("Yield (a.u.)", fontsize=16)
ax1.set_xlabel("Harmonic order", fontsize=16)
ax1.tick_params(which="major",axis='both',labelsize=16)
ax1.set_title("HHG spectrum", fontsize=20)
ax1.set_ylim(1E-7,1e5)

fig.tight_layout()
plt.savefig(f"{path}/Spectrum.png")
plt.show()
plt.close()
