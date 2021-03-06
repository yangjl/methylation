#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import argparse
import textwrap
import timeit
import os

def version():
    ver0 = """
    ##########################################################################################
    splitcontext version 0.1.0
    Author: Jinliang Yang
    purpose: Split the data by context
    --------------------------------
    
    updated: 09-02-2016, first piece of the code
    ##########################################################################################
    """
    return ver0


    

##########################################################################################
def read_write(infile, outbase):
    
    cgfile = ".".join([outbase, "cg"])
    chgfile = ".".join([outbase, "chg"])
    chhfile = ".".join([outbase, "chh"])
    cg, chg, chh, cg1, chg1, chh1 = (0,0,0,0,0,0)
    
    with open(infile, 'r') as infile, open(cgfile, 'w') as outcg, open(chgfile, 'w') as outchg, open(chhfile, 'w') as outchh:

        for line in infile:
            token = line.split()
            if(token[5] == "CG"):
                cg += 1
                if((int(token[3]) + int(token[4])) != 0):
                    cg1 += 1
                    outcg.write("\t".join(token[0:5]) + "\n")
            elif(token[5] == "CHG"):
                chg += 1
                if((int(token[3]) + int(token[4])) != 0):
                    chg1 += 1
                    outchg.write("\t".join(token[0:5]) + "\n")
            elif(token[5] == "CHH"):
                chh += 1
                if((int(token[3]) + int(token[4])) != 0):
                    chh1 += 1
                    outchh.write("\t".join(token[0:5]) + "\n")
            else:
                print("WARNING: !!! unknow context !!!", token[5])
                sys.exit()
    
    ### output statistics
    outstat = ".".join([outbase, "stat"])
    with open(outstat, 'w') as outs:
        outs.write("\t".join(["context", "methylated", "total"]) + "\n")
        outs.write("\t".join(["CG", str(cg1), str(cg)]) + "\n")
        outs.write("\t".join(["CHG", str(chg1), str(chg)]) + "\n")
        outs.write("\t".join(["CHH", str(chh1), str(chh)]) + "\n")

#f='/Users/yangjl/Desktop/simusnps.txt'
#read_write(infile="largedata/test_CX_report.txt", outbase="largedata/test_res")
 
####    
def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )

    # positional arguments:
    #parser.add_argument('query', metavar='QUERY', type=str, nargs='*', \
    #                    help='the question to answer')

    # optional arguments:
    parser.add_argument('-p', '--path', help='the path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-i','--input', help='input file: generated by bismark', type=str)
    parser.add_argument('-o', '--output', help='output files', type=str)

    return parser
    #parser = get_parser()
    #parser.print_help()

if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['input'] is not None:
        print(version())
    if args['path'] is not None:
        os.chdir(args['path'])

    ##### cal running time ######
    st = timeit.default_timer()
    
    print(">>> Reading and writing data ...")
    ### checking the input file
    read_write(infile=args['input'], outbase=args['output'])

    et = timeit.default_timer()

    print(">>> [ ", "%.0f" % ((et - st)/60), " ] minutes of run time!")
    print(">>> Job finished!")
