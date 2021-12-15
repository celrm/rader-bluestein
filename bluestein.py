import math
import random

from fft_pot import fft,ifft

def bluestein(a):
    #Calculos basicos
    n = len(a)
    xi = pow(math.e,1j * math.tau/n)

    #Calculamos u(n) y v(n)
    u = [a[i] * pow(xi,-(i*i)/2) for i in range (n)]
    v = [pow(xi,(i *i)/2) for i in range (n)]
    v_estrella = [pow(xi,-(i * i)/2) for i in range (n)]

    #Extension a potencia de 2
    l = 2 ** math.ceil(math.log2(2 * n + 1))
    xi_l = math.cos(math.tau / l) + 1j * math.sin(math.tau / l)

    #Extension de u
    u_l = u + [0] * (l - n)
    

    #Extension de v
    aux = v[1:]
    aux.reverse()
    v_l = v + [0] * (l - 2*n + 1) + aux

    #Calculo de la convolucion

    dft_u_l = fft(u_l, xi_l)  
    dft_v_l = fft(v_l, xi_l)  

    uv_l = [i*j for i,j in zip(dft_u_l,dft_v_l)]

    ift_l = ifft(uv_l , xi_l) 

    ift = ift_l[:n]

    # multiplicamos por el factor v*
    dft = [i*j for i,j in zip(v_estrella,ift)]

    return dft


# Comparar con la versión O(n²)
def check_dft(a, dft):
    dft_correcta = [
        sum([a[i] * pow(math.e, -math.tau * 1j * j * i / len(a)) for i in range(len(a))])
        for j in range(len(a))
    ]
    print(sum([abs(dft[i] - dft_correcta[i]) for i in range(len(a))]))


p = 11
a = [random.uniform(0, 1) for _ in range(p)]
check_dft(a, bluestein(a))
