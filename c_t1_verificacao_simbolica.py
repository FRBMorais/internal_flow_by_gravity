from sympy import symbols, log, simplify

# Definindo as variáveis
r, k, mu, g_z, nu, r_i, c1v = symbols('r k mu g_z nu r_i c1v')

# Expressões simbólicas anteriores
c1t_simb = mu*(16*c1v**2*nu**2*log(r_i) + 8*c1v*g_z*nu*r_i**2 + g_z**2*r_i**4)/(16*k*nu**2)

# Sua expressão para c1t
c1t_num = (mu/k) * ((g_z ** 2 * r_i ** 4 / (16 * nu ** 2)) + (g_z * c1v * r_i ** 2 / (2 * nu)) + c1v ** 2 * log(r_i))

# Diferença entre os termos c1t
delta_c1t = simplify(c1t_simb - c1t_num)

print("Diferença simplificada entre os termos c1t:")
print(delta_c1t)

# okokokokokokoko