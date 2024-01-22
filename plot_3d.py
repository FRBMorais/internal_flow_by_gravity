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

for nome, prop in FLUIDOS.items():
    nu = prop['nu']
    mu = prop['mu']
    k = prop['k']
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    axs = axs.flatten()

    for idx, dif in enumerate(diferencas):
        re = ri + dif
        N = 100
        dr = (re - ri) / (N-1)
        r = np.linspace(ri, re, N)
        
        c1v = -gz * (re**2 - ri**2) / (4 * nu * np.log(re / ri))
        
        A = np.zeros((N, N))
        b = np.zeros(N)

        for i in range(1, N - 1):
            rp = r[i] + dr/2
            rm = r[i] - dr/2

            A[i, i - 1] = k * rp / (dr**2 * r[i])
            A[i, i] = -k * (rp + rm) / (dr**2 * r[i])
            A[i, i + 1] = k * rm / (dr**2 * r[i])

            dvz_dr = gz * r[i] / (2 * nu) + c1v / r[i]
            b[i] = -mu * dvz_dr**2

        # Condição de contorno de Neumann na parede interna
        A[0, 0] = 1
        A[0, 1] = -1

        # Condição de contorno de convecção na parede externa
        A[N - 1, N - 1] = -k / dr - h
        A[N - 1, N - 2] = k / dr
        b[N - 1] = -h * t_infinito

        T = np.linalg.solve(A, b)

        axs[idx].plot(r, T, label=f'Raio Externo: {re} m')
        axs[idx].set_title(f'Diferença: {dif} m')
        axs[idx].set_xlabel('Raio (m)')
        axs[idx].set_ylabel('Temperatura (C)')
        axs[idx].legend()

    fig.suptitle(nome)
    plt.show()
