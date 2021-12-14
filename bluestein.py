
def bluestein(a,xi):
    #Calculamos u(n) y v(n)
    n = len(a)
    u = [a[i] * pow(xi, pow(i,2)/2) for i in n]
    v = [pow(-xi, pow(i,2)/2) for i in n]
    v_estrella = [pow(xi, pow(i,2)/2) for i in n]

    #DFT u
    d_u = fft(u,xi)
    #DFT v
    d_v = fft(v,xi)

    uv = [i*j for i,j in zip(d_u,d_v)]

    i_uv = ifft(uv,xi)

    dft = [i*j for i,j in zip(v_estrella,i_uv)]

    return dft
