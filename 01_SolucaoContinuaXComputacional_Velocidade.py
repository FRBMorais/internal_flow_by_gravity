import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

G_Z = 9.81


def equacao_diferencial(r, Vz, propriedades):
    """Equação diferencial para resolver usando BVP.

    Parâmetros:
        r (np.ndarray): Array representando o raio.
        Vz (np.ndarray): Array representando a velocidade.
        propriedades (dict): Propriedades do fluido.

    Retorna:
        list: Valores calculados com base na equação diferencial.
    """
    nu = propriedades['nu']
    dVz = Vz[1]
    d2Vz = G_Z / nu - (1 / r) * dVz
    return [dVz, d2Vz]


def condicoes_de_contorno(Vz_a, Vz_b):
    """Condições de contorno para a equação diferencial.

    Parâmetros:
        Vz_a (np.ndarray): Condições iniciais de velocidade.
        Vz_b (np.ndarray): Condições finais de velocidade.

    Retorna:
        list: Lista das condições de contorno.
    """
    return [Vz_a[0], Vz_b[0]]


def perfil_de_velocidade_analitico(r, nu, r_i, r_e):
    """Solução analítica para o perfil de velocidade.

    Parâmetros:
        r (np.ndarray): Array representando o raio.
        nu (float): Viscosidade cinemática.
        r_i (float): Raio inicial.
        r_e (float): Raio externo.

    Retorna:
        np.ndarray: Perfil de velocidade para os valores de raio fornecidos.
    """
    c1v = (G_Z * (r_i ** 2 - r_e ** 2)) / (4 * nu * np.log(r_e / r_i))
    c2v = (-G_Z / (4 * nu * np.log(r_e / r_i))) * (r_i ** 2 * np.log(r_e / r_i) + (r_i ** 2 - r_e ** 2) * np.log(r_i))
    return (G_Z * r ** 2 / (4 * nu)) + c1v * np.log(r) + c2v, c1v


def plotar_graficos(fluidos):
    """Gera e exibe gráficos para os fluidos fornecidos.

    Parâmetros:
        fluidos (dict): Dicionário de propriedades dos fluidos.
    """
    r_i = 0.05
    diferencas = [0.005, 0.010, 0.015, 0.020]
    conjuntos_de_raios = [(r_i, r_i + diff) for diff in diferencas]

    for nome_fluido, propriedades in fluidos.items():
        fig_velocidade, axs_velocidade = plt.subplots(2, 2, figsize=(15, 10))
        axs_velocidade = axs_velocidade.ravel()
        fig_erro, axs_erro = plt.subplots(2, 2, figsize=(15, 10))
        axs_erro = axs_erro.ravel()

        for idx, (r_i, r_e) in enumerate(conjuntos_de_raios):
            r = np.linspace(r_i + 1e-6, r_e - 1e-6, 1000)
            Vz_init = np.zeros((2, r.size))

            # Solução numérica
            sol = solve_bvp(lambda r, Vz: equacao_diferencial(r, Vz, propriedades), condicoes_de_contorno, r, Vz_init)
            Vz_numerico = sol.sol(r)[0]

            # Solução analítica
            Vz_analitico, _ = perfil_de_velocidade_analitico(r, propriedades['nu'], r_i, r_e)

            # Cálculos de erro
            erro_absoluto = np.abs(Vz_numerico - Vz_analitico)
            erro_relativo = erro_absoluto / (np.abs(Vz_analitico) + 1e-10)

            # Plot do perfil de velocidade
            axs_velocidade[idx].plot(r, Vz_numerico, 'b-', label='Numérico')
            axs_velocidade[idx].plot(r, Vz_analitico, 'r--', label='Analítico')
            axs_velocidade[idx].set_xlabel('Raio r (m)')
            axs_velocidade[idx].set_ylabel('Velocidade Vz (m/s)')
            axs_velocidade[idx].set_title(f'Comparação com $r_e - r_i = {r_e - r_i:.2f}$ m')
            axs_velocidade[idx].legend()
            axs_velocidade[idx].grid(True)

            # Plot de erro
            axs_erro[idx].plot(r, erro_absoluto, 'g-.', label='Erro Absoluto')
            axs_erro[idx].plot(r, erro_relativo, 'm-.', label='Erro Relativo')
            axs_erro[idx].set_xlabel('Raio r (m)')
            axs_erro[idx].set_ylabel('Erro')
            axs_erro[idx].set_title(f'Erros com $r_e - r_i = {r_e - r_i:.2f}$ m')
            axs_erro[idx].legend()
            axs_erro[idx].grid(True)

        fig_velocidade.suptitle(f'Comparação de Soluções para {nome_fluido}')
        fig_velocidade.tight_layout()
        fig_velocidade.subplots_adjust(top=0.88)

        fig_erro.suptitle(f'Erros para {nome_fluido}')
        fig_erro.tight_layout()
        fig_erro.subplots_adjust(top=0.88)

    plt.show(block=True)


if __name__ == "__main__":
    FLUIDOS = {
        'Água': {'nu': 0.0000008007, 'mu': 0.0007972, 'k': 0.598},  # ok
        'Ar': {'nu': 0.00001589, 'mu': 0.00001846, 'k': 0.0263},  # ok
        'Óleo de motor': {'nu': 0.00055, 'mu': 0.486, 'k': 0.145},  # ok
        'Glicerina': {'nu': 0.000634, 'mu': 0.799, 'k': 0.286}  # ok
    }

    plotar_graficos(FLUIDOS)
