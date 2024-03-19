import numpy as np

#  a1 and a3 are, respectively, the intensity and wavenumber,obtained in the quantum chemical calculations.
# The remaining shape parameters, a2 and a4, were set as 0.075 and 0.015, following the best agreement with the experimental bandshapes.
def A(v, a3 , a1, a2=0.075, a4=0.015):
  return a1 * np.exp(-(a4**2) * (v - a3)**2) / (1 + a2**2 * (v - a3)**2)

data = {
    1573.34829012: 76.49504446,
    3109.58706778: 0.59058089,
    3648.74297346: 2.04505505,
    3730.39770829: 2.46557077,
    5234.41957248: 0.10609848,
    5279.99719162: 0.00484495,
    7224.01921962: 0.01027934,
    7240.01703666: 0.04340794,
    7380.32774933: 0.00570194,
}

# for wavenumber, intensity in data.items():
#   result = A(0, wavenumber, intensity )
#   print(result)

# exit()

ls = [0]*10001
for wavenumber, intensity in data.items():
  for i in range(10000):
    result = A(i,wavenumber,intensity)
    ls[i] = ls[i] + result
    if abs(i-wavenumber) < 1 :
      print(i, wavenumber, intensity, result, ls[i])


print(ls)

import matplotlib.pyplot as plt

plt.plot(ls)
plt.xlabel('Index')
plt.ylabel('Intensity')
plt.show()
