import sys, os, re

def reformat(f):
    of = 'reformatted/'+f
    with open(f, 'r') as input, open(of, 'w') as output:
        line = input.readline() # Comment 1
        output.write(line)
        line = input.readline() # Comment 2
        output.write(line)
        line = input.readline() # Comment 3
        output.write(line)

        line = input.readline() # Number of elements and elements
        output.write(line)
        sline = line.strip().split()
        while '' in sline:
            sline.remove('')
        nelements = int(sline[0])
        
        line = input.readline() # nrho drho nr dr cutoff
        output.write(line)
        sline = line.strip().split()
        while '' in sline:
            sline.remove('')
        nrho = int(sline[0])
        drho = float(sline[1])
        nr = int(sline[2])
        dr = float(sline[3])
        cutoff = float(sline[4])

        line = 'holder'
        while line != []:
            line = input.readline().strip().split()
            while '' in line:
                line.remove('')
            try:
                [float(x) for x in line]
                for num in line:
                    output.write(num+'\n')
            except:
                output.write(' '.join(line) + '\n')


def main():
    f = sys.argv[1]
    if not os.path.exists(os.path.split(os.path.abspath(f))[0] + '/reformatted'):
        os.makedirs(os.path.split(os.path.abspath(f))[0] + '/reformatted')
    reformat(f)

if __name__ == '__main__':
    main()
