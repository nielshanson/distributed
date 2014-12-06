#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2010, The metapaths Project"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     from os import makedirs, sys, remove, path
     import re
     from optparse import OptionParser, OptionGroup

     from LCAComputation import *
     from metapaths_utils  import  fprintf, printf
     #from libs.python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)


usage= """./tax.py --tax-file taxfile """

parser = OptionParser(usage)
parser.add_option("--tax-file", dest="tax_file",
                  help='name of the refseq names file')


def copyList(a, b):
    [ b.append(x) for x in a ]


def get_species(hit):
    species = []
    try:
        m = re.findall(r'\[([^\[]+)\]', hit)
        if m != None:
          copyList(m,species)
    except:
          return None

    if species:
       return species
    else:
       return None


def create_annotation(namefile):

    file = 'ncbi_taxonomy_tree.txt'
    lca = LCAComputation(file)

    #taxonomy=lca.getTaxonomy(species)


def read_taxonomy_list(filename,namesdict):
    namefile = open(filename, 'r')
    lines = namefile.readlines()
  
    for line in lines:
       line = line.strip()
       if line:
          names = get_species(line)
          if names:
            for name in names:
                namesdict[name]=1

   
# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    namesdict = {}

   
    read_taxonomy_list(opts.tax_file, namesdict)

    file = 'ncbi_taxonomy_tree.txt'
    lca = LCAComputation(file)
    lca.construct_tree()
    sys.exit(0)

    for name in namesdict:
       print '=' + name
       lca.get_lineage(name)
       

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

