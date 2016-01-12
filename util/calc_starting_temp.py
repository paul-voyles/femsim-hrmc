import sys, math, random
import numpy as np

rand = lambda: 1.0-random.random()

def main():
    # Import a file containing values for del-chi for the first ~1000 moves
    acceptance = 0.90
    kb = 8.6171e-05
    jobfile = sys.argv[1]

    content = []
    start = False
    i = 0
    for line in open(jobfile):
        if 'Starting step' in line:
            start = True
        if start:
            if 'Del-chi' in line:
                i += 1
                x = float(line.strip().split('=')[1])
                #if x < 0: continue
                #if x == 0: x = 1e-14
                content.append(x)
            if i > 50000:
                break
 
    #delfile = sys.argv[1]
    #content = open(delfile).readlines()
    #content = [float(line.strip()) for line in content]

    #values = [math.exp(-x/(kb*T)) for x in content]
    values = [x for x in content]
    T = 50
    dt = 10

    I = 100
    J = 300
    i = 0
    j = 0
    of = open('temperature.txt', 'w')
    avgs = [0 for i in range(J)]
    while True:
        rs = [rand() for x in values]
        compare = [-kb*T*math.log(r) for r in rs]
        val = sum(compare[i] > values[i] for i in range(len(values)))
        val = float(val)/len(values)
        avgs[j] = T #val/(1+abs(acceptance-val))
        print(T, val, j, sum(avgs)/len(avgs), dt)
        of.write("{0}  {1}\n".format(T, val))
        if val < acceptance:
            T += dt
        elif val > acceptance:
            T -= dt
        else:
            break
        if i >= I-1:
            dt = dt/2.0
            i = -1
        if j >= J-1:
            j = -1
        i += 1
        j += 1
    print(T)
    of.close()


if __name__ == '__main__':
    main()
