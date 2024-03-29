def finndata(filnavn="log.lammps",l0="Step Time",l1="Loop time",verbose=False):
    with open(filnavn,"r") as fil:
        linjer = fil.readlines()
    N = len(linjer)
    data = {}

    for i in range(N):
        l = linjer[i]
        if l0 in l:
            if verbose:
                print(l)
            titler = l.split()
            for ord in titler:
                if not ord in data.keys():
                    data[ord] = []
            i += 1
            while i<N and not l1 in linjer[i]:
                ord = linjer[i].split()
                for j in range(len(ord)):
                    data[titler[j]].append(float(ord[j]))
                i += 1
    return data

def smooth(y,bredde):
    # Bredde = antall paa EN side
    from numpy import asarray,zeros,mean
    y = asarray(y)
    N = len(y)
    ny = zeros(N)
    for i in range(bredde,N-bredde):
        ny[i] = mean(y[i-bredde:i+bredde+1])
    for i in range(bredde):
        ny[i] = mean(y[:bredde])
        ny[-i-1] = mean(y[N-bredde:N])
    return ny


if __name__=="__main__":
    from matplotlib.pyplot import *
    from argparse import *

    parser = ArgumentParser("Skript som plotter filer med log.lammps-format.")

    parser.add_argument("-x","--x",default="Step",dest="x")
    parser.add_argument("-y","--y",default="c_C_mass_center[3]",dest="y")
    parser.add_argument("-f","--fil",nargs="+",dest="f")
    parser.add_argument("-l0","--startdatalinje",default="Step Time",dest="l0")
    parser.add_argument("-l1","--sluttdatalinje",default="Loop time",dest="l1")
    parser.add_argument("-s","--save",default=None,dest="s")
    parser.add_argument("--smooth",default=0,type=int)
    parser.add_argument("--dump",default=None)
    parser.add_argument("--dumpnum",default=0,type=int)

    args = parser.parse_args()

    x = []
    y = []

    if not type(args.f) == list:
        args.f = [args.f]

    for fil in args.f:
        data = finndata(fil,args.l0,args.l1,verbose=True)
        x += data[args.x]
        y += data[args.y]

    if not args.smooth == 0:
        y = smooth(y,args.smooth)


    finished_length = min(len(x),len(y))
    plot(x[:finished_length],y[:finished_length])
    xlabel(args.x)
    ylabel(args.y)
    grid()
    if not args.s is None:
        savefig(args.s,bbox_inches="tight")
    show()
    if not args.dump is None:
        with open(args.dump,"w") as fil:
            N=len(x)
            if args.dumpnum == 0:
                args.dumpnum = N
            for i in range(0,N,N/args.dumpnum):
                fil.write("%g %g\n" % (x[i],y[i]))
