# De clase: FFT para potencias de 2
def fft(pol, xi):
    n = len(pol)
    if n == 1:
        return pol
    pol_par = [pol[2 * i] for i in range(n // 2)]
    pol_impar = [pol[2 * i + 1] for i in range(n // 2)]
    res_par = fft(pol_par, xi ** 2)
    res_impar = fft(pol_impar, xi ** 2)
    res = [0j] * n
    for i in range(n // 2):
        res[i] = res_par[i] + xi ** i * res_impar[i]
        res[i + n // 2] = res_par[i] - xi ** i * res_impar[i]
    return res


# De clase: Inversa de FFT para potencias de 2
def ifft(vec, xi):
    n = len(vec)
    res = fft(vec, 1 / xi)
    for i in range(n):
        res[i] /= n
    return res
