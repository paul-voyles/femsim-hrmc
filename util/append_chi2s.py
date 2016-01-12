import sys,os
import natsort

def main():
    path = sys.argv[1]
    files = os.listdir(path)
    if('chi2_total.txt' in files):
        files.remove('chi2_total.txt')
    files = [f for f in files if 'chi' in f]
    files = natsort.natsorted(files)

    of = open('chi2_total.txt','w')
    of.write('step, chi2, energy\n')
    for f in files:
        inf = open(f)
        inf.readline()
        for line in inf:
            of.write(line)
        inf.close()
    of.close()


if __name__ == "__main__":
    main()
