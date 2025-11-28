import matplotlib.pyplot as plt
import numpy as np

#
x = np.linspace(-0.1, 0.1, 200)

y1 = 0.8 + 0.03*np.sin(20*x)
y2 = np.exp(-(x/0.03)**2) * 0.45
y3 = np.exp(-(x/0.05)**2) * 0.35
y4 = np.exp(-(x/0.06)**2) * 0.20
y5 = np.full_like(x, 0.35)

plt.figure(figsize=(5.5, 5.5))

plt.plot(x, y1, color="black", linewidth=2)
plt.plot(x, y2, color="magenta", linewidth=2)
plt.plot(x, y3, color="orange", linewidth=2)
plt.plot(x, y4, color="brown", linewidth=2)
plt.plot(x, y5, color="grey", linewidth=2)

plt.ylim(0, 1)
plt.title("S1")

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_visible(True)

plt.xticks([-0.1, 0, 0.1], ["-0.1", "CEN", "+0.1"])
plt.grid(False)
plt.xlabel("")
plt.ylabel("")

plt.tight_layout()
plt.show()