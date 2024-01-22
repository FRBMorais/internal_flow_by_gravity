import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Definindo variáveis
r, nu, g_z, r_i, r_e, mu, rho = sp.symbols('r nu g_z r_i r_e mu rho')
c_1v, c_2v = sp.symbols('c_1v c_2v')

# Expressão para V_z(r)
V_z = g_z * r**2 / (4*nu) + c_1v * sp.ln(r) + c_2v

# Aplicando condições de contorno
boundary_cond_1 = V_z.subs(r, r_i)
boundary_cond_2 = V_z.subs(r, r_e)

# Resolvendo para c_1v e c_2v usando as condições de contorno
sol = sp.solve([boundary_cond_1, boundary_cond_2], (c_1v, c_2v))

# Valores fornecidos para c_1v e c_2v
c1v_given = -g_z * (r_e**2 - r_i**2) / (4*nu*sp.ln(r_e/r_i))
c2v_given = -g_z * (r_i**2 * sp.ln(r_e) - r_e**2 * sp.ln(r_i)) / (4*nu*sp.ln(r_e/r_i))

# Simplificando a diferença
diff_c1v_simplified = sp.simplify(sol[c_1v] - c1v_given)
diff_c2v_simplified = sp.simplify(sol[c_2v] - c2v_given)

# Substituindo os valores de c_1v e c_2v encontrados na expressão V_z
V_z = V_z.subs({c_1v: sol[c_1v], c_2v: sol[c_2v]})

# Calculando dV_z/dr
dVz_dr = sp.diff(V_z, r)

# Potência de transformação viscosa
P_v = 2 * sp.pi * sp.integrate(mu * r * dVz_dr**2, (r, r_i, r_e))

# Velocidade média dada a igualdade entre Pv e Pg
V = P_v / (rho * g_z)

# Calculando a potência mecânica gravitacional usando V
P_g = rho * g_z * V

# Verificando a igualdade
diff_power = sp.simplify(P_g - P_v)

# Calculando a velocidade média integrando o campo de velocidade
integral_Vz = 2 * sp.pi * sp.integrate(V_z * r, (r, r_i, r_e))
area_total = sp.pi * (r_e**2 - r_i**2)
V_integral = integral_Vz / area_total

# Verificar se as velocidades médias são iguais
diff_velocities = sp.simplify(V - V_integral)

print(f"Potencia de transformacao viscosa == \n{P_v}\n"
      f"Potencia mecanica gravitacional == \n{P_g}\n")

if diff_power == 0:
    print("A potência mecânica gravitacional e a potência de transformação viscosa são iguais.\n")
else:
    print("As potências não são iguais.\n")
    print(f"Diferença entre P_g e P_v: \n{diff_power}\n")

print(f"Velocidade média calculada pela relação Pv = rho g_z V: \n{V}\n")
print(f"Velocidade média calculada pela integral do campo de velocidade: \n{V_integral}\n")

if diff_velocities == 0:
    print("As duas velocidades médias são iguais.\n")
else:
    print(f"Diferença entre as duas velocidades médias: \n{diff_velocities}\n\n")

# Imprimindo os resultados para constantes
print(f"c_1v_computacional --> \n{sol[c_1v]}\n")
print(f"c_1v_manual --> \n{c1v_given}\n")
print(f"Diferença simplificada entre c_1v manual e calculado: \n{diff_c1v_simplified}")

print(f"\nc_2v_computacional --> \n{sol[c_2v]}\n")
print(f"c_2v_manual --> \n{c2v_given}\n")
print(f"Diferença simplificada entre c_2v manual e calculado: \n{diff_c2v_simplified}")

