import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parâmetros fornecidos
delta_ri = 0.01
r_inicial = 1.5
r_e = 2.5
vel_ang_i = 144
vel_ang_e = 60

# Parâmetros da discretização
r_values = np.arange(r_inicial, r_e + delta_ri, delta_ri)
n_r = len(r_values)
dt = 0.001
t_final = 10  # tempo total de simulação

# Inicialização do campo de velocidade
V_theta = np.zeros(n_r)
V_theta[0] = vel_ang_i
V_theta[-1] = vel_ang_e


# Função para atualizar V_theta usando método explícito
def update_V_theta(V_theta, dt, delta_ri):
    new_V = V_theta.copy()
    for k in range(1, n_r - 1):
        r_k = r_inicial + k * delta_ri
        r_k_half = r_inicial + (k + 0.5) * delta_ri
        r_k_minus_half = r_inicial + (k - 0.5) * delta_ri

        diffusive_term = (r_k_half * V_theta[k + 1] - (r_k_half + r_k_minus_half) * V_theta[k] + r_k_minus_half *
                          V_theta[k - 1]) / delta_ri ** 2
        angular_term = -V_theta[k] / r_k

        new_V[k] = V_theta[k] + dt * (diffusive_term + angular_term)

    return new_V


# Animação
fig, ax = plt.subplots()
line, = ax.plot(r_values, V_theta, label="V_theta")
ax.set_ylim(0, 150)
ax.set_xlabel('r')
ax.set_ylabel('V_theta')
ax.legend()


def animate(i):
    global V_theta
    V_theta = update_V_theta(V_theta, dt, delta_ri)
    line.set_ydata(V_theta)
    ax.set_title(f"Time: {i * dt:.2f}s")
    return line,


ani = animation.FuncAnimation(fig, animate, frames=int(t_final / dt), blit=True, repeat=False)
plt.show()
