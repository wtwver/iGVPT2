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

with open("PET.out", "r") as f:
    lines = f.readlines()
    result = False
    result_ls = []
    for idx, line in enumerate(lines, 1):
        if "Results With intensities" in line:
            result = not result
        if result:
            line = line.strip().split(" ")
            line = [elem for elem in line if elem.strip()]
            result_ls.append(line)


result_dict = {float(i[1]): float(i[3]) for i in result_ls[5:-3]}


ls = [0]*10001
for wavenumber, intensity in result_dict.items():
  for i in range(10000):
    result = A(i,wavenumber,intensity)
    ls[i] = ls[i] + result
    #if abs(i-wavenumber) < 1 :
      #print(i, wavenumber, intensity, result, ls[i])


import matplotlib.pyplot as plt

plt.plot(ls)
plt.xlabel('Index')
plt.ylabel('Intensity')
plt.show()