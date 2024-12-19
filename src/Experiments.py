import pandas as pd
from matplotlib import pyplot as plt

df = pd.read_excel("Test.xlsx")
data = df.to_numpy()[:, 0:8]
print(data[0, :])


# REMAINDER REGIONS AS A FUNCTION OF N
N = data[:-1, 1]
t = data[:-1, 2]
d = data[:-1, 6]
Rr = data[:-1, 4]
h = data[:-1, 3]
n = data[:-1, 0]

RrT = 2 * N - 7 + 7 * h

plt.scatter(N, Rr, color='blue', label='Number of Remainder Regions (r)', s=5)
plt.scatter(N, RrT, color='red', label='Theoritical Bound', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('r (y-axis)')
plt.title('The number of Remainder Regions (r) as a function of N')
plt.legend()
plt.show()

plt.scatter(N, RrT-Rr, color='blue',
            label='The difference between r and the theoretical bound', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('r-theoretical bound (y-axis)')
plt.title('The difference between r and the theoretical bound as a function of N')
plt.legend()
plt.show()


# DISKS AS A FUNCTION OF N
dT = 2 * N - 7 + 7 * h
plt.scatter(N, d, color='blue', label='Number of Disks (d)', s=5)
plt.scatter(N, dT, color='red', label='Theoritical Bound', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('d (y-axis)')
plt.title('The number of Disks (d) as a function of N')
plt.legend()
plt.show()

plt.scatter(N, dT-d, color='blue',
            label='The difference between d and the theoretical bound', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('d-theoretical bound (y-axis)')
plt.title('The difference between d and the theoretical bound as a function of N')
plt.legend()
plt.show()

plt.scatter(n, t/n, color='blue',
            label='The number of Triangles (t) divided by n', s=5)
plt.xlabel('n (x-axis)')
plt.ylabel('t/n (y-axis)')
plt.title('The number of Triangles (t) dividided by n as a function of n')
plt.legend()
plt.show()

plt.scatter(t/Rr, Rr/RrT, color='blue',
            label='The number of Triangles per Remainder Region', s=5)
plt.xlabel('The ratio of the number of Remainder Regions (x-axis)')
plt.ylabel('The number of Triangles per Remainder Region (y-axis)')
plt.title('The impact of the ratio of the number of Remainder Regions')
plt.legend()
plt.show()

plt.scatter(N, t/N, color='blue',
            label='The number of Triangles (t) divided by N', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/N (y-axis)')
plt.title('The number of Triangles (t) dividided by N as a function of N')
plt.legend()
plt.show()

plt.scatter(N, t/Rr, color='blue',
            label='The number of Triangles divided by Remainer Regions (r)', s=5)
plt.xlabel('N (x-axis)')
plt.ylabel('t/Rr (y-axis)')
plt.title(
    'The number of Triangles divided by Remainer Regions (r) as a function of N')
plt.legend()
plt.show()
