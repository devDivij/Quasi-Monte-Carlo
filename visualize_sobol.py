import matplotlib.pyplot as plt

x = []
y = []

with open("sobol.txt") as f:
    for line in f:
        a,b = map(float, line.split())
        x.append(a)
        y.append(b)

plt.figure(figsize=(5,5))
plt.scatter(x, y, s=10)
plt.title("Sobol Points in 2D")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()
