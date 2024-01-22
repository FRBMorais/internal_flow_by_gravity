import numpy as np
import matplotlib.pyplot as plt

FLUIDOS = {
    'Água': {'nu': 1.004e-6, 'mu': 0.001, 'k': 0.606},
    'Ar': {'nu': 1.48e-5, 'mu': 1.85e-5, 'k': 0.0257},
    'Óleo de motor': {'nu': 1e-4, 'mu': 0.09, 'k': 0.145},
    'Glicerina': {'nu': 0.92e-2, 'mu': 1.49, 'k': 0.286}
}

ri = 0.5
diferencas = [0.01, 0.05, 0.1, 0.5]
gz = 9.81
t_infinito = 30
h = 20

class TemperatureField:
    def __init__(self, g_z, nu, micro, k, h, T_infinito):
        self.g_z = g_z
        self.nu = nu
        self.micro = micro
        self.k = k
        self.h = h
        self.T_infinito = T_infinito

    def c_1v(self, r_i, r_e):  # ok
        return -self.g_z * (r_e**2 - r_i**2) / (4 * self.nu * np.log(r_e / r_i))

    def c_t1(self, r_i, c1v):  # ok
        return (self.micro / self.k) * (self.g_z**2 * r_i**4 / (16 * self.nu**2) + c1v**2 * np.log(r_i) + self.g_z * c1v * r_i**2 / (2 * self.nu))

    def c_t2(self, r_i, r_e, c1v, c1t):  # ok
        term1 = self.T_infinito - c1t * np.log(r_e) - c1t * self.k / (r_e * self.h)
        term2 = (self.micro / self.h) * (self.g_z * r_e**3 / (16 * self.nu**2) + self.g_z * c1v * r_e / (2 * self.nu) + c1v**2 * np.log(r_e) / r_e)
        term3 = (self.micro / self.k) * (self.g_z * r_e**4 / (64 * self.nu**2) + self.g_z * c1v * r_e**2 / (4 * self.nu) + c1v**2 * np.log(r_e)**2 / 2)
        return term1 + term2 + term3

    def t_r(self, r, c1v, c1t, c2t):  # ok
        term1 = (-self.micro / self.k) * (self.g_z * r**4 / (64 * self.nu**2) + self.g_z * c1v * r**2 / (4 * self.nu) + c1v**2 * np.log(r)**2 / 2)
        return term1 + c1t * np.log(r) + c2t


for nome, prop in FLUIDOS.items():
    nu = prop['nu']
    mu = prop['mu']
    k = prop['k']
    
    campo_temperatura = TemperatureField(gz, nu, mu, k, h, t_infinito)
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    axs = axs.flatten()

    for idx, dif in enumerate(diferencas):
        re = ri + dif
        N = 100
        dr = (re - ri) / (N-1)
        r = np.linspace(ri, re, N)
        
        A = np.zeros((N, N))
        b = np.zeros(N)

        for i in range(1, N - 1):
            rp = r[i] + dr / 2
            rm = r[i] - dr / 2
            A[i, i - 1] = k * rp / (dr**2 * r[i])
            A[i, i] = -k * (rp + rm) / (dr**2 * r[i])
            A[i, i + 1] = k * rm / (dr**2 * r[i])

            dvz_dr = gz * r[i] / (2 * nu) + campo_temperatura.c_1v(ri, re) / r[i]
            b[i] = -mu * dvz_dr**2

        A[0, 0] = 1
        A[0, 1] = -1
        A[N - 1, N - 1] = -k / dr - h
        A[N - 1, N - 2] = k / dr
        b[N - 1] = -h * t_infinito

        T = np.linalg.solve(A, b)
        
        c1v = campo_temperatura.c_1v(ri, re)
        c1t = campo_temperatura.c_t1(ri, c1v)
        c2t = campo_temperatura.c_t2(ri, re, c1v, c1t)
        T_analitico = [campo_temperatura.t_r(raio, c1v, c1t, c2t) for raio in r]

        axs[idx].plot(r, T, label=f'Raio Externo: {re} m - Numérico')
        axs[idx].plot(r, T_analitico, '--', label=f'Raio Externo: {re} m - Analítico')
        axs[idx].set_title(f'Diferença: {dif} m')
        axs[idx].set_xlabel('Raio (m)')
        axs[idx].set_ylabel('Temperatura (C)')
        axs[idx].legend()

    fig.suptitle(nome)
    plt.show()
