import math
from convolution import conv
import random

def bluestein(a):
    n = len(a)
    xi = math.cos(math.tau / n) + 1j * math.sin(math.tau / n)
    #Calculamos u(n) y v(n)
    
    u = [a[i] * pow(xi, pow(i,2)/2) for i in range(n)]
    v = [pow(xi, -pow(i,2)/2) for i in range(n)]
    v_estrella = [pow(xi, pow(i,2)/2) for i in range(n)]

    u_conv_v = conv(u, v)

    dft = [i*j for i,j in zip(v_estrella,u_conv_v)]

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
check_dft(a, bluestein(a))