import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Dimensões do cilindro em corte
altura_cilindro = 5.0
largura_cilindro = 1.0

# Espessura da camada de líquido fino
espessura_liquido = 0.1

fig, ax = plt.subplots()

# Desenhando o cilindro em corte
cilindro = patches.Rectangle((0, 0), largura_cilindro, altura_cilindro, linewidth=1, edgecolor='k', facecolor='none', label='Cilindro')
ax.add_patch(cilindro)

# Desenhando a camada líquida fino nas bordas laterais
liquido_fino_esquerda = patches.Rectangle((0, 0), espessura_liquido, altura_cilindro, linewidth=1, edgecolor='k', facecolor='blue', alpha=0.3)
ax.add_patch(liquido_fino_esquerda)

liquido_fino_direita = patches.Rectangle((largura_cilindro - espessura_liquido, 0), espessura_liquido, altura_cilindro, linewidth=1, edgecolor='k', facecolor='blue', alpha=0.3)
ax.add_patch(liquido_fino_direita)

# Vetores de raios
ax.annotate('', xy=(espessura_liquido, altura_cilindro/2), xytext=(largura_cilindro/2, altura_cilindro/2), arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=1.5), va='center')
ax.text(espessura_liquido/2 + espessura_liquido, altura_cilindro/2 - 0.2, r'$r_i$', verticalalignment='center', horizontalalignment='center', fontsize=12)

ax.annotate('', xy=(largura_cilindro, altura_cilindro/2), xytext=(largura_cilindro/2, altura_cilindro/2), arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=1.5), va='center')
ax.text(largura_cilindro - espessura_liquido / 2, altura_cilindro/2 - 0.2, r'$r_e$', verticalalignment='center', horizontalalignment='center', fontsize=12)

ax.annotate('', xy=(largura_cilindro, altura_cilindro/2), xytext=(largura_cilindro/2, altura_cilindro/2), arrowprops=dict(facecolor='black', arrowstyle='-|>', lw=1.5), va='center')
ax.text(largura_cilindro - espessura_liquido / 2, altura_cilindro/2 - 0.2, r'$r$', verticalalignment='center', horizontalalignment='center', fontsize=12)

# Propriedades do fluido
props_fluido = r'Propriedades do fluido interno'
ax.text(largura_cilindro + 0.1, altura_cilindro/2 + 0.8, props_fluido, verticalalignment='center', horizontalalignment='left', fontsize=12, color='red')

props_fluido = r'$\mu$, $k$, $\rho$, $c_p$'
ax.text(largura_cilindro + 0.5, altura_cilindro/2 + 0.5, props_fluido, verticalalignment='center', horizontalalignment='left', fontsize=12, color='red')

# Propriedades do ar externo
props_ar_externo = r'Propriedades do fluido externo'
ax.text(largura_cilindro + 0.1, altura_cilindro/2 - 0.2, props_ar_externo, verticalalignment='center', horizontalalignment='left', fontsize=12, color='red')

# Propriedades do ar externo
props_ar_externo = r'$h$, $T_{\infty}$'
ax.text(largura_cilindro + 0.5, altura_cilindro/2 - 0.5, props_ar_externo, verticalalignment='center', horizontalalignment='left', fontsize=12, color='red')

# Opções de visualização
ax.set_xlim([-1, largura_cilindro + 2])
ax.set_ylim([0, altura_cilindro + 1])
ax.set_title("Corte Longitudinal do Cilindro")
plt.show()
