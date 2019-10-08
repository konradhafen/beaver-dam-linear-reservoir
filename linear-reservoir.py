import numpy as np
import matplotlib.pyplot as plt


def getArea(h0, h1, s, p):
    area = s * (h1/h0) ** (2/p)
    return area

def getHeight(h0, v, s, p):
    h = (h0**(2*p)*((v+(2*v)/p)/s))**(1/(1+2*p))
    return h

def getVolume(h0, h, s, p):
    v = (s/(1+(2/p)))*((h**(1+(2/p))/(h0**(2*p))))
    return v


v0 = 1135.0  # initial volume (m^3)
a0 = 2974.0  # initial area (m^2)
s = 889.0  # area scaling factor (m^2)
h0 = 1.0  # maximum height/pond depth (m)
p = 0.69  # represents shape of pond's bathymetric curve (1=linear, >1=concave, <1=convex)
etr = 0.00  # et rate (m/day)
k = 0.2  # reservoir/recession coefficient (day^-1)

nponds = 2
ndays = 5
volume = np.zeros((ndays, nponds))
inflow = np.copy(volume)
outflow = np.copy(volume)
et = np.copy(volume)
area = np.copy(volume)

volume[0, :] = v0
inflow[:, 0] = np.linspace(1.0*24*60*60/35.3147, 0, ndays)
h = np.full((1, nponds), h0)  # inital height (water depth) for each pond

for i in range(0, volume.shape[0]-1):
    area[i, :] = getArea(h0, h, s, p)
    et[i, :] = etr * area[i, :]
    outflow[i, :] = volume[i, :] * k
    volume[i, :] = volume[i, :] - et[i, :] - outflow[i, :]
    inflow[i, 1:] = inflow[i, 1:] + outflow[i, :-1]
    excess = v0 - (volume[i, :] + inflow[i, :] - et[i, :] - outflow[i, :])
    volume[i+1, :] = volume[i, :] + inflow[i, :] - et[i, :] - outflow[i, :]
    excess = volume[i + 1, :] - v0
    print("outflow 1", outflow[i, :])
    outflow[i, :] = outflow[i, :] + np.where(excess > 0, excess, 0.0)
    print("outflow 2", outflow[i, :])
    volume[i + 1, :][excess > 0] = v0
    inflow[i, 1:] = inflow[i, 1:] + outflow[i, :-1]
    volume[volume < 0] = 0.0
    h = getHeight(h0, volume[i, :], s, p)
outflow[-1, :] = volume[-1, :] * k


plt.plot(np.linspace(0, ndays, ndays), (outflow*35.3147)/(24*60*60))
plt.show()
plt.plot(np.linspace(0, ndays, ndays), np.sum(et, axis=1))
plt.show()
