import matplotlib.pyplot as plt

# Listas para armazenar os valores de X e Y de cada arquivo
allX = []
allY = []

# Cores para cada plotagem
colors = ["red", "green", "blue", "orange", "purple", "magenta", "cyan", "yellow"]

# Nome da instancia

instance = "B"

# Nomes dos arquivos
paths = [f"{instance}/elisama.txt", f"{instance}/original.txt", f"{instance}/islame.txt", f"{instance}/original2.txt", f"{instance}/islame2.txt", f"{instance}/elisama3.txt", f"{instance}/original3.txt", f"{instance}/islame3.txt"]

# Função para ler os arquivos e preencher as listas allX e allY
def read_file(file_name):
    X, Y = [], []
    with open(file_name, 'r') as f:
        for line in f:
            if len(line.split()) == 0:
                # Quando encontrar uma linha vazia, salva os valores e limpa as listas temporárias
                allX.append(X.copy())
                allY.append(Y.copy())
                X.clear()
                Y.clear()
            else:
                # Adiciona os valores de cada linha
                x, y = line.split()
                X.append(float(x))
                Y.append(float(y))
        # Adiciona o último conjunto de valores ao final da leitura
        if X and Y:
            allX.append(X.copy())
            allY.append(Y.copy())

# Lê os três arquivos
for file in paths:
    read_file(file)

plt.title(f"{instance}")
plt.xlabel("Cost (max)")
plt.ylabel("Power (max)")

# Plota os dados de cada conjunto de pontos
for i in range(len(paths)):
    plt.plot(allX[i], allY[i], color='none', label=paths[i], marker='o', markersize=4, markerfacecolor=colors[i])
 
# Mostra a legenda e o gráfico
plt.legend()
plt.show()
