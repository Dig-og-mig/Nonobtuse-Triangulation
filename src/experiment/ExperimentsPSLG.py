import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

df = pd.read_excel("dataPSLG.xlsx")
data = df.to_numpy()[:, 0:8]
data = np.delete(data, 6, axis=0)
print(data[:, :])

df2 = pd.read_excel("Test.xlsx")
data2 = df2.to_numpy()[:, 0:8]
Np = data2[:-1, 1]
tp = data2[:-1, 2]
Rrp = data2[:-1, 4]


# REMAINDER REGIONS AS A FUNCTION OF N
N = data[:, 1]
t = data[:, 2]
d = data[:, 6]
Rr = data[:, 4]
h = data[:, 3]
n = data[:, 0]
n2 = n**2.5

# plt.scatter(n, t/n, color='blue', label='The number of Triangles (t) divided by n', s = 5)
# plt.xlabel('n (x-axis)')
# plt.ylabel('t/n (y-axis)')
# plt.title('The number of Triangles (t) dividided by n as a function of n')
# plt.legend()
# plt.show()


plt.scatter(n, t/n, color='blue',
            label='The number of Triangles (t) divided by n', s=5)
plt.ylim(0, 800)
plt.xlabel('n (x-axis)')
plt.ylabel('t/n (y-axis)')
plt.title('The number of Triangles (t) dividided by n as a function of n (PSLG)')
plt.legend()
plt.show()

plt.scatter(n, t, color='blue', label='The number of Triangles (t)', s=5)
# plt.ylim(0, 800)
plt.xlabel('n (x-axis)')
plt.ylabel('t (y-axis)')
plt.title('The number of Triangles (t) as a function of n (PSLG)')
plt.legend()
plt.show()

plt.scatter(n, t/n2, color='blue',
            label='The number of Triangles (t) divided by n^2.5', s=5)
plt.ylim(0, 50)
plt.xlabel('n (x-axis)')
plt.ylabel('t/n^2.5 (y-axis)')
plt.title('The number of Triangles (t) dividided by n^2.5 as a function of n (PSLG)')
plt.legend()
plt.show()
plt.scatter(n, t/n2, color='blue',
            label='The number of Triangles (t) divided by n^2.5', s=5)
plt.ylim(0, 12.5)
plt.xlabel('n (x-axis)')
plt.ylabel('t/n^2.5 (y-axis)')
plt.title('The number of Triangles (t) dividided by n^2.5 as a function of n (PSLG)')
plt.legend()
plt.show()

plt.scatter(N, t/N, color='blue',
            label='The number of Triangles (t) divided by N', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/N (y-axis)')
plt.title('The number of Triangles (t) dividided by N as a function of N (PSLG)')
plt.legend()
plt.show()

plt.scatter(Np, tp/Np, color='blue',
            label='The number of Triangles (t) divided by N', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/N (y-axis)')
plt.title('The number of Triangles (t) dividided by N as a function of N (polygon)')
plt.legend()
plt.show()

plt.scatter(N, t/N, color='red',
            label='The number of Triangles (t) divided by N (PSLG)', s=5)
plt.scatter(Np, tp/Np, color='blue',
            label='The number of Triangles (t) divided by N (polygon)', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/N (y-axis)')
plt.title(
    'The number of Triangles (t) dividided by N as a function of N (polygon + pslg)')
plt.legend()
plt.show()

plt.scatter(N, t/Rr, color='red',
            label='The number of Triangles (t) divided by r (PSLG)', s=5)
plt.scatter(Np, tp/Rrp, color='blue',
            label='The number of Triangles (t) divided by r (polygon)', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/N (y-axis)')
plt.title(
    'The number of Triangles (t) dividided by r as a function of N (polygon + pslg)')
plt.legend()
plt.show()
