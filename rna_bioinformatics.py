from math import factorial
import nussinov
import numpy as np
import matplotlib.pyplot as plt


# attachment 1

def maximum_matching_rna(seq: str) -> float:
    a = seq.count("A")
    u = seq.count("U")
    g = seq.count("G")
    c = seq.count("C")

    au = factorial(max(a, u)) / factorial(max(a, u) - min(a, u))
    gc = factorial(max(g, c)) / factorial(max(g, c) - min(g, c))

    return au * gc


print(maximum_matching_rna("CCGG"))  # returns 2
print(maximum_matching_rna("CCGGG"))  # returns 6
print(maximum_matching_rna("CCGGAAUUU"))  # returns 12

# attachment 2

k = 0.001985875  # kcal/mol*K
T = 273 + 37  # Kelvin
en = [-28.10, -27.90, -27.80, -27.80, -27.60, -27.50, -27.20, -27.20, -27.20, -27.20, -27.10, -27.00, -27.00, -27.00,
      -27.00, -27.00, -27.00, -26.70, -26.70, -26.70, -26.70, -26.70]
x = sorted(list(set(en)))
Z = 0
y = []
y2 = []
print(len(list(x)))
d = dict(zip(x, [en.count(el) for el in x]))
print(d)
for i in en:
    y.append(np.exp(-i / (k * T)))
for i in x:
    y2.append(d[i] * np.exp(-i / (k * T)))
Z = sum(y)
y = np.array(y) / Z
y2 = np.array(y2) / Z
plt.plot(en, y)
plt.show()

x_en = list()
y_en = list()
for energy, y in d.items():
    x_en.append(energy)
    y_en.append(y * np.exp(-energy / (k * T)))

energy_y_sum = sum(y_en)
y_en = np.array(y_en) / energy_y_sum

summed = plt.bar(x_en, y_en, align="center", fill=False)
plt.show()


