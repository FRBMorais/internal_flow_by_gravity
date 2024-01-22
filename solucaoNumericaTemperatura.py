import numpy as np
import matplotlib.pyplot as plt

# Constantes
g_z = 9.81  # m/s^2
nu = 1e-6  # m^2/s
r_i = 0.5  # m

# Criando 100 valores de r_e
r_es = np.linspace(0.5005, 0.55, 100)

plt.figure(figsize=(10, 6))

# Loop sobre os diferentes valores de r_e
for r_e in r_es:
    c_1v = (g_z * (r_i ** 2 - r_e ** 2)) / (4 * nu * np.log(r_e / r_i))
    
    c_2v = (-g_z / (4 * nu * np.log(r_e / r_i))) * (r_i ** 2 * np.log(r_e / r_i) + (r_i ** 2 - r_e ** 2) * np.log(r_i))

    # Criando um array de valores de r entre r_i e r_e
    r = np.linspace(r_i, r_e, 100)

    # Calculando a derivada da velocidade para cada valor de r
    print(f'parcela linear: {(g_z * r / (2 * nu))}')
    print(f'parcela hiperbolica: {c_1v / r}')
    velocity_gradients = (g_z * r / (2 * nu)) + c_1v / r

    # Plotando o gráfico
    plt.plot(r, velocity_gradients)

plt.xlabel('Raio $r$ (m)')
plt.ylabel('Derivada da Velocidade $\\frac{dV}{dr}$ $(s^{-1})$')
plt.title('Derivada da Velocidade para diferentes valores de $r_e$')
plt.grid(True)
# A linha abaixo pode ser descomentada se você quiser ver as legendas, mas com 100 curvas pode ficar poluído.
# plt.legend()
plt.show()
