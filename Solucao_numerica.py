import numpy as np
import matplotlib.pyplot as plt

# Constantes
g_z = 9.81  # m/s^2, aceleração da gravidade
nu = 1e-6  # m^2/s, viscosidade cinemática da água
micro = 1e-3  # Pa·s, viscosidade dinâmica da água
k = 0.6  # W/(m*K), condutividade térmica da água
r_i = 0.5  # m, raio interno

# Criando 100 valores de r_e
r_es = np.linspace(0.5005, 0.55, 100)

plt.figure(figsize=(10, 6))

# Loop sobre os diferentes valores de r_e
for r_e in r_es:
    # Calculando os coeficientes c_1v e c_1t para cada valor de r_e
    c_1v = (g_z * (r_i ** 2 - r_e ** 2)) / (4 * nu * np.log(r_e / r_i))
    c_1t = (micro / k) * (
                g_z ** 2 * r_i ** 4 / (16 * nu ** 2) + g_z * c_1v * r_i ** 2 / (2 * nu) + c_1v ** 2 * np.log(r_i))

    # Criando um array de valores de r entre r_i e r_e
    r = np.linspace(r_i, r_e, 100)

    # Calculando o gradiente de temperatura para cada valor de r
    gradients = -micro / k * (
                g_z ** 2 * r ** 3 / (16 * nu ** 2) + g_z * c_1v * r / (2 * nu) + c_1v ** 2 * np.log(r) / r) + c_1t / r

    # Plotando o gráfico
    plt.plot(r, gradients, label=f"$r_e$ = {r_e:.4f}")

plt.xlabel('Raio $r$ (m)')
plt.ylabel('$\\frac{dT}{dr}$ (K/m)')
plt.title('Gradiente de Temperatura para diferentes valores de $r_e$')
plt.grid(True)
# A linha abaixo pode ser descomentada se você quiser ver as legendas, mas com 100 curvas pode ficar poluído.
# plt.legend()
plt.show()
