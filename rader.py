import math
import random
from convolution import conv


# Potencias sucesivas: [1, gen, gen², ...]
def secuencia(gen, p):
    return [pow(gen, i, p) for i in range(p - 1)]


# Generador de (Z/pZ)* y sus potencias sucesivas (p primo)
def generador(p):
    for gen in range(2, p):
        perm = secuencia(gen, p)
        if len(perm) == p - 1:
            return gen, perm


# MÉTODO Rader: FFT de a
def rader(a):
    p = len(a)  # p primo

    # g, [1, g, g²,...]
    gen, gen_seq = generador(p)
    # g⁻¹, [1, g⁻¹, g⁻²,...]
    inv = pow(gen, p - 2, p)
    inv_seq = secuencia(inv, p)

    # u(n) = ξ^g⁻ⁿ
    xi_p = math.cos(math.tau / p) + 1j * math.sin(math.tau / p)  # ξ raíz p-ésima
    u = [pow(xi_p, g_n) for g_n in inv_seq]
    # v(m) = a_gᵐ
    v = [a[g_m] for g_m in gen_seq]
    u_conv_v = conv(u, v)

    dft = [0] * p
    dft[0] = sum(a)  # dft₀ = Σᵢ aᵢ
    for j in range(p - 1):  # dft_g⁻ʲ = a₀ + (u * v)ⱼ
        dft[inv_seq[j]] = a[0] + u_conv_v[j]
    return dft


# Comparar con la versión O(n²)
def check_dft(a, dft):
    dft_correcta = [
        sum([a[i] * pow(math.e, math.tau * 1j * j * i / len(a)) for i in range(len(a))])
        for j in range(len(a))
    ]
    print(sum([abs(dft[i] - dft_correcta[i]) for i in range(len(a))]))


p = 11
a = [random.uniform(0, 1) for _ in range(p)]
check_dft(a, rader(a))
