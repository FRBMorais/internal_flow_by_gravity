import numpy as np
import matplotlib.pyplot as plt


class TemperatureField:
    def __init__(self, g_z, nu, micro, k, h, T_infinito):
        self.g_z = g_z
        self.nu = nu
        self.micro = micro
        self.k = k
        self.h = h
        self.T_infinito = T_infinito

    def c_1v(self, r_i, r_e):  # ok
        return -self.g_z * (r_e ** 2 - r_i ** 2) / (4 * self.nu * np.log(r_e / r_i))  # ok

    def c_t1(self, r_i, c1v):  # ok
        return (self.micro / self.k) * (
                self.g_z ** 2 * r_i ** 4 / (16 * self.nu ** 2) + c1v ** 2 * np.log(r_i) + self.g_z * c1v * r_i ** 2 / (
                2 * self.nu))

    def c_t2(self, r_i, r_e, c1v, c1t):
        term1 = self.T_infinito - c1t * np.log(r_e) - c1t * self.k / (r_e * self.h)  # ok

        term2 = (self.micro / self.h) * (
                self.g_z ** 2 * r_e ** 3 / (16 * self.nu ** 2) + self.g_z * c1v * r_e / (2 * self.nu) + c1v ** 2 * np.log(
            r_e) / r_e)  # ok

        term3 = (self.micro / self.k) * (
                self.g_z ** 2 * r_e ** 4 / (64 * self.nu ** 2) + self.g_z * c1v * r_e ** 2 / (4 * self.nu) + c1v ** 2 * np.log(
            r_e) ** 2 / 2)  # ok

        return term1 + term2 + term3

    def t_r(self, r, c1v, c1t, c2t):
        term1 = (-self.micro / self.k) * (
                self.g_z ** 2 * r ** 4 / (64 * self.nu ** 2) + self.g_z * c1v * r ** 2 / (4 * self.nu) + c1v ** 2 * np.log(
            r) ** 2 / 2)
        return term1 + c1t * np.log(r) + c2t


FLUIDOS = {
    'Água': {'nu': 1.004e-6, 'mu': 0.001, 'k': 0.606},
    'Ar': {'nu': 1.48e-5, 'mu': 1.85e-5, 'k': 0.0257},
    'Óleo de motor': {'nu': 1e-4, 'mu': 0.09, 'k': 0.145},
    'Glicerina': {'nu': 0.92e-2, 'mu': 1.49, 'k': 0.286}
}

r_i = 0.5
diferencas = [0.01, 0.05, 0.1, 0.5]
g_z = 9.81
t_infinito = 300
h = 20


def equacao(r, re, fluido):
    nu = FLUIDOS[fluido]['nu']
    k = FLUIDOS[fluido]['k']
    c_1v = -g_z * (re ** 2 - r_i ** 2) / (4 * nu * np.log(re / r_i))

    dVz_dr = g_z * r / (2 * nu) + c_1v / r
    micro = FLUIDOS[fluido]['mu']

    d2Tdr2 = -(k / micro) * dVz_dr ** 2

    return d2Tdr2


def resolver_e_plotar(mostrar_numerico=True, mostrar_continuo=True):
    for fluido in FLUIDOS.keys():
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle(f'Temperatura em função de r para {fluido}')

        axs = axs.flatten()

        for i, dif in enumerate(diferencas):
            r_e = r_i + dif
            r = np.linspace(r_i, r_e, 100)
            T = np.zeros_like(r)

            temp_field = TemperatureField(g_z, FLUIDOS[fluido]['nu'], FLUIDOS[fluido]['mu'], FLUIDOS[fluido]['k'], h, t_infinito)
            c1v = temp_field.c_1v(r_i, r_e)
            c1t = temp_field.c_t1(r_i, c1v)
            c2t = temp_field.c_t2(r_i, r_e, c1v, c1t)
            T_cont = [temp_field.t_r(ri, c1v, c1t, c2t) for ri in r]

            for j in range(1, len(r) - 1):
                dr = r[j] - r[j - 1]
                d2Tdr2 = equacao(r[j], r_e, fluido)
                T[j + 1] = 2 * dr ** 2 * d2Tdr2 + 2 * T[j] - T[j - 1]

            k = FLUIDOS[fluido]['k']
            T[-1] = (h * t_infinito + k * T[-2]) / (h + k)

            if mostrar_numerico:
                axs[i].plot(r, T, label=f'Numérica - r_e = {r_e}')
            if mostrar_continuo:
                axs[i].plot(r, T_cont, label=f'Contínua - r_e = {r_e}', linestyle='dashed')

            axs[i].set_title(f'r_e = {r_e}')
            axs[i].set_xlabel('r')
            axs[i].set_ylabel('Temperatura (°C)')
            axs[i].legend()
            axs[i].grid(True)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()

# Exemplo de uso:
resolver_e_plotar(mostrar_numerico=True, mostrar_continuo=False)

