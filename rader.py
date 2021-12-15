import math
import random
from fft_pot import fft, ifft  # funciones hechas en clase


# Potencias sucesivas: [1, gen, gen², ...]
def secuencia(gen, p):
    return [pow(gen, i, p) for i in range(p - 1)]


# Generador de (Z/pZ)* y sus potencias sucesivas (p primo)
def generador(p):
    for gen in range(2, p):
        perm = secuencia(gen, p)
        if len(perm) == p - 1:
            return gen, perm


# Convolución de u y v
# Relleno de 0s hasta l (potencia de 2)
def conv(u, v):
    n = len(u)  # l = 2ᵏ > 2n
    l = 2 ** math.ceil(math.log2(2 * n + 1))
    xi_l = math.cos(math.tau / l) + 1j * math.sin(math.tau / l)  # ξ raíz l-ésima

    dft_ul = fft(u + [0] * (l - n), xi_l)  # fft de u extendido
    dft_vl = fft(v + [0] * (l - n), xi_l)  # fft de v extendido

    # (uℓ * vℓ) = ifft(dft(uℓ) · dft(vℓ))
    ul_conv_vl = ifft([a * b for (a, b) in zip(dft_ul, dft_vl)], xi_l)

    # reducir a tamaño n de nuevo
    return [ul_conv_vl[j] + ul_conv_vl[j + n] for j in range(n)]


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
