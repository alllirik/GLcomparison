import sys
import getopt
from matplotlib import pyplot as plt
import pandas as pd

def myfunc(argv):
    arg_coverage = ""
    arg_help = "{0} -c <coverage>".format(argv[0])
    
    try:
        opts, args = getopt.getopt(argv[1:], "hc:", ["help", "coverage="])
    except:
        print(arg_help)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-c", "--coverage"):
            arg_coverage = arg

    cov=arg_coverage
    print('coverage:', arg_coverage)
    sfs=[] #This list will contain sfs estimates from different models
    for name in ['truth.sfs','samtools.sfs','gatk.sfs','atlas.sfs']:
        f = open(cov+'.'+name)
        res=''
        for i in f:
            res+=i
        res=[float(x) for x in res.strip().split()]
        print(res)
        sfs.append(res)
        
    #trueSFS * chrLen     
    #for i in range(len(sfs[0])):
    #    sfs[0][i]*=1000000

    n=len(sfs[0])

    sfs_df=pd.DataFrame({'ATLAS':sfs[0][1:n-1], 'SAMtools':sfs[1][1:n-1], 'GATK':sfs[2][1:n-1], 'ATLAS':sfs[3][1:n-1]})
    sfs_df.index = sfs_df.index + 1 #Allele frequency stars with 1
    fig = plt.figure(figsize=(10, 6), dpi=90)
    ax = fig.add_subplot(1,1,1) 
    width = 0.15

    sfs_df.Truth.plot(kind='bar', color='black', ax=ax, width=width, position=0)
    sfs_df.SAMtools.plot(kind='bar', color='red', ax=ax, width=width, position=-1)
    sfs_df.GATK.plot(kind='bar', color='green', ax=ax, width=width, position=-2)
    sfs_df.ATLAS.plot(kind='bar', color='blue', ax=ax, width=width, position=-3)

    ax.set_ylabel('Number of sites')
    ax.set_xlabel('Derived allele frequency')
    ax.legend(prop={'size': 12})
    plt.title('Site Frequency Spectrum: ' + str(int((n-1)/2)) + ' Simulated individuals')
    plt.savefig('simulated_'+cov+'_SFS.png')

if __name__ == "__main__":
    myfunc(sys.argv)

