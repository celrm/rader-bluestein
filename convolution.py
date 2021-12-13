import math
from fft_pot import fft, ifft  # funciones hechas en clase

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
