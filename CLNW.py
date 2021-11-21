import argparse
from random import randint
import time
from sage.all import *

A = 87450056743442738320404100044938503254010774441639884906892991297288383737801
B = 58300037828961825546936066696625668836007182961093256604595327531525589158534
p = 115792089237316193816632940749697632406671544878850700875012406725593640484581
x = 52452215485807069130338229023131158406785328425500997187572698503722848469538
y = 34390556094101336083038232990243039188087353598375496864910039276079318512317
k = 333333333300000000003333333333
# log_2(k) ~ 98: Optimal window length = 4

# 111111111110111111111111111111111111111111101111111111111111111011111111111111111111111111111111111111111111111111111110111111111    - 680398580342365126779515082699991678463
# 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001    - 340282366920938463463374607431768211457

# Calculation method: CLNW
def clnw(E, k, P, pr):
    Q = E([0, 1, 0])

    k_b = bin(k)[2:]
    t = int(ceil(log(k, 2)))
    d = windowsSize(t)
    
    if (d == 0):
        print("\nWrong windows size.")
        exit(-1)
    
    if (pr == 1):
        print(f'\nWindows size = {d}')
        print(f'k = {k_b}')
    
    stop = pow(2, d)
    points = []

    ptime = time.time()
    points.append(0)
    for i in range(1, stop):
        if((i % 2) == 1):
            points.append(P * i)
        else:
            points.append(0)
    ptime = time.time() - ptime

    split_k = windows(k, d)
    s = len(split_k)
    int_sk = []
    len_sk = []
    
    if (pr == 1):
        print('Windows:', end = ' ')
    for i in range(0, s, 1):
            int_sk.append(int(split_k[i], 2))
            len_sk.append(len(split_k[i]))
            if (pr == 1):
                if (int_sk[i] != 0):
                    print(f'NZW_{i} = ({split_k[i]})', end = ' ')
                else:
                    print(f'ZW_{i} = ({split_k[i]})', end = ' ')
    if (pr == 1):
        print('\n')

    rtime = time.time()
    for i in range(s - 1, -1, -1):
        Q *= pow(2, len_sk[i])
        j = int_sk[i]
        if (j != 0):
            Q += points[j]
    rtime = time.time() - rtime
    
    if (pr == 1):
        sagePk = P * k
        print(f'k * P = {k} * {P} = {Q}')
        print(f'Sage: k * P = {k} * {P} = {sagePk}')
        
        return (ptime, rtime, (sagePk == Q))

    return (ptime, rtime, Q)

def windowsSize(t):
    if (t < 128): 
        return 4
    if (t >= 128 and t < 768):
        return 5
    if (t >= 768 and t < 1792):
        return 6
    if (t >= 1792):
        return 7

    return 0

def windows(k, d):
    k_s = []
    zw = ''
    nzw = ''
    t = bin(k)[2:]
    j = len(t) - 1
    
    while(j >= 0):
        if (t[j] == '0'):
            if(len(nzw) > 0):
                k_s.append(nzw)
                nzw = ''
            zw += '0'
            j -= 1
        else:
            if(len(zw) > 0):
                k_s.append(zw)
                zw = ''
            if(len(nzw) > 0):
                k_s.append(nzw)
                nzw = ''
            if (j >= d):
                for i in range (j - d + 1, j + 1, 1):
                   nzw += t[i]
                j = j - d
            else:
                for i in range(j, -1, -1):
                   nzw += t[i]
                j = j - d
                i = len(nzw)
                for y in range(i, d):
                   nzw = '0' + nzw
                k_s.append(nzw)
        
    return k_s

# Binary algorithm for calculating a multiple point
def kPBinCalculate(E, k, P, pr):
    Q = E([0, 1, 0])
    k_t = bin(k)[2:]

    t = int(ceil(log(k, 2)))

    rtime = time.time()
    for i in range(t):
        Q = 2 * Q
        if (k_t[i] == '1'):
            Q = Q + P    
    rtime = time.time() - rtime

    if (pr == 1):
        sagePk = P * k
        print(f'\nk * P = {k} * {P} = {Q}')
        print(f'Sage: k * P = {k} * {P} = {sagePk}')
        
        return (rtime, (sagePk == Q))

    return (rtime, Q)

def experiments(E, i, k):
    bin_time = 0
    clnw_time = 0
    
    while (i > 0):
        i -= 1
        P = E.random_point()

        rtime, Q = kPBinCalculate(E, k, P, 0)
        # print(f'k * P = {k} * {P} = {Q}')
        bin_time += rtime
    
        ptime, rtime, Q = clnw(E, k, P, 0)
        # print(f'k * P = {k} * {P} = {Q}')
        clnw_time += rtime
    
    return (bin_time, clnw_time)
    
    
def start(E, xp, yp, k, operation):    
    try:
        P = E([xp, yp])
    except:
        print(u'ERROR: Wrong point P coordinates. A random point on the curve is selected.')
        return (-2)
    else:
        if (operation == 'calcCLNW' or operation == 'calcV3' or operation == 'clnw'):
            p_time, r_time, res = clnw(E, k, P, 1)

            if res:
                print(f'\nRuntime for preprocessing multiplications = {p_time} seconds')
                print(f'Runtime for algorithm = {r_time} seconds')
                print(f'Total runtime for algorithm = {r_time + p_time} seconds\n')
            else:
                print("Wrong result.")
                return (-1)

        elif (operation == 'calculate' or operation == 'calcBin' or operation == 'clc'):
            r_time, res = kPBinCalculate(E, k, P, 1)

            if res:
                print(f'\nRuntime = {r_time} seconds\n')
            else:
                print("Wrong result.")
                return (-1)

        else:
            print (u'ERROR: Need to select an operation.')
            return (-1)

    return (1)

def main():
    Info = argparse.ArgumentParser()
    Info.add_argument("-p", help = 'field characteristic', type = int, default = p)
    Info.add_argument("-A", help = 'coefficient A of the curve equation in Weierstrass form', type = int, default = A)
    Info.add_argument("-B", help = 'coefficient B of the curve equation in Weierstrass form', type = int, default = B)
    Info.add_argument("-x", help = 'x-coordinate of point P', type = int, default = x)
    Info.add_argument("-y", help = 'y-coordinate of point P', type = int, default = y)
    Info.add_argument("-k", help = 'number k', type = int, default = k)
    Info.add_argument("-i", help = 'number of iterations for the experiments', type = int, default = 10)
    Info.add_argument("-o", help = 'operation - [calcCLNW, calcV3, clnw] / [calculate, calcBin, clc] / [experiments, exp]')
    
    InfoParsed = Info.parse_args()
    new_p = InfoParsed.p
    new_A = InfoParsed.A
    new_B = InfoParsed.B
    xp = InfoParsed.x
    yp = InfoParsed.y
    new_k = InfoParsed.k
    i = InfoParsed.i
    operation = InfoParsed.o
        
    try:
        E = EllipticCurve(GF(new_p), [new_A, new_B])
    except:
        print(u'ERROR: Wrong elliptic curve parameters.')
        Info.print_help()
        exit(-1)
    else:
        if (operation == 'experimentse' or operation == 'exp'):
            for j in range (1, 2):
                b_time, c_time = experiments(E, i**j, new_k)
            
                print(f'\nNumber of iterations = {i**j}')
                print(f'Runtime of binary algorithm = {b_time} seconds')
                print(f'Runtime of CLNW algorithm  = {c_time} seconds\n')
        
        else:
            st = start(E, xp, yp,  new_k, operation)
            if (st == -2):
                P = E.random_point()
                st = start(E, P[0], P[1], new_k, operation)
            if (st == -1):
                Info.print_help()
                exit(-1)

if __name__ == '__main__':
    main()