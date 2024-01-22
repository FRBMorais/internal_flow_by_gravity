import numpy as np
import matplotlib.pyplot as plt

fluids = {
    'Água': {'nu': 1.004e-6, 'mu': 0.001, 'k': 0.606},
    'Ar': {'nu': 1.48e-5, 'mu': 1.85e-5, 'k': 0.0257},
    'Óleo de motor': {'nu': 1e-4, 'mu': 0.09, 'k': 0.145},
    'Glicerina': {'nu': 0.92e-2, 'mu': 1.49, 'k': 0.286}
}

# Variáveis
r_i_val = 0.05
r_e_values = np.linspace(r_i_val, 0.07, 20)[1:]  # Começando do segundo valor
delta_r_values = r_e_values - r_i_val  # Diferença entre r_e e r_i
g_z_val = 9.81  # Aceleração da gravidade


# Função para calcular a potência de transformação viscosa
def potencia_transformacao_viscosa(mu, g_z, nu, r_e, r_i):
    return mu * g_z ** 2 * np.pi / (8 * nu ** 2) * (
                (r_e ** 4 - r_i ** 4) - ((r_e ** 2 - r_i ** 2) ** 2 / np.log(r_e / r_i)))


# Plot
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.ravel()

for idx, (fluid_name, properties) in enumerate(fluids.items()):
    P_v_values = [potencia_transformacao_viscosa(properties['mu'], g_z_val, properties['nu'], r_e_val, r_i_val) for
                  r_e_val in r_e_values]

    axs[idx].plot(delta_r_values, P_v_values, label=f"P_v {fluid_name}")
    axs[idx].set_title(fluid_name)
    axs[idx].set_xlabel('delta r |r_e - r_i| [m]')
    axs[idx].set_ylabel('Potência de Transformação Viscosa [W/m]')
    axs[idx].legend(loc="upper left")
    axs[idx].text(0.1, 0.1, f'r_i = {r_i_val} [m]', transform=axs[idx].transAxes,
                  bbox=dict(facecolor='white', alpha=0.5))
    axs[idx].grid(True)

plt.tight_layout()
plt.show()
