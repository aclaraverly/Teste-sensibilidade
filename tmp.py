import pandas as pd
import matplotlib.pyplot as plt

# Carregar o CSV usando delimitador ";" e decimal ","
df = pd.read_csv("article_data.csv", delimiter=";", decimal=",")

# Renomear para clareza (opcional)
df.columns = ["x", "y"]

# Fazer o gráfico: eixo x = segunda coluna, eixo y = primeira
plt.figure(figsize=(8,6))
plt.plot(df["x"], df["y"], marker="o", linestyle="-")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Gráfico de y em função de x")
plt.grid(True, linestyle="--", alpha=0.7)
plt.show()