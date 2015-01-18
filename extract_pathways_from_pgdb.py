#!/usr/bin/python

# import libraries
import sys # for argument vector
import optparse # to parse arguments
import pythoncyc # pythoncyc
import operator # sorting hashes
import pickle
import os

# describe what the script does
what_i_do = "A simple script to access pathway tools"
usage = "extract_pathways_from_pgdb.py --ptools_organism meta -o output_file"

# initialize the parser
parser = optparse.OptionParser(usage = usage, epilog = what_i_do)
parser.add_option("-p", "--ptools_organism", dest="ptools_org", default=None,
                   help='Organism to attach pathway tools in [Required]')
parser.add_option("-l", "--list_organisms", action="store_true", dest="list_organisms", default=None,
                  help='List available PGDBs in Pathway Tools.')                   
parser.add_option("-o", "--output_file", dest="output_file", default=None,
                   help='output file to write to [Required]')


# check required options
def check_arguments(opts):
    if not opts.ptools_org:
        print "Error: Pathway Tools organism/pgdb id not specified."
        print usage
        exit()
    # if not opts.output_file:
    #         print "Error: Output file not specified."
    #         print usage
    #         exit()

# connect to pathway tools organism fail otherwise
def connect_to_ptools(ptools_org):
    try:
        cyc = pythoncyc.select_organism(ptools_org)
    except:
        print "Error: Could not connect to ptools"
        exit()
    return cyc

# cleans off list 
def clean_ptools_output(my_list):
    new_list = []
    for x in my_list: 
        new_list.append(x.strip("|"))
    return new_list

# Tests for edge between two reactions. If reactions share an edge then they should
# be connected 
def test_for_edge(rxn1, rxn2):
    # collect reactants and product
    l1 = set(cyc.get_slot_values(rxn1, "LEFT"))
    r1 = set(cyc.get_slot_values(rxn1, "RIGHT"))
    l2 = set(cyc.get_slot_values(rxn2, "LEFT"))
    r2 = set(cyc.get_slot_values(rxn2, "RIGHT"))
    
    substrate_list = list(l1) + list(r1) + list(l2) + list(r2)
    for item in substrate_list:
        if item not in metabolite_count:
            metabolite_count[item] = 0
        metabolite_count[item] += 1
    
    if (l1.intersection(r2) or r1.intersection(l2)):
        # intersection not empty
        return True
    else:
        return False

# adds an edge to the graph
def add_edge(source, target, rxn, pwy=None, edge_type="path"):
    # could be a problem
    new_edge = "_".join([source, target, rxn, pwy, edge_type])
    if new_edge not in edges:
        

# add edges to graph from connnected component from pathway pwy
def add_edges_from_component(component, pwy):
    for rxn in component[::-1]:
           left = cyc[rxn]['left']
           right = cyc[rxn]['right']
           for l in left:
               for r in right:
                   if (l not in excluded_compounds) and (r not in excluded_compounds):
                       add_edge(l, r, rxn, pwy, "component")
    for pwy in all_edges:
        print pwy
        for rxn in all_edges[pwy]:
            print "\t", rxn
            for source in 

# the main function
def main():
   (opts, args) = parser.parse_args()
   check_arguments(opts)
   
   # connect to pathway tools
   global cyc
   cyc = pythoncyc.select_organism(opts.ptools_org)
   
   global excluded_compounds # list of uninformative compounds
   excluded_compounds = {}
   
   if os.path.exists("/tmp/metabolite_count.pk"):
       fh = open("/tmp/metabolite_count.pk", 'r')
       metabolite_count = pickle.load(fh)
       fh.close()
       sorted_metabolite_count = sorted(metabolite_count.items(), key=operator.itemgetter(1), reverse=True)
       for metabolite in sorted_metabolite_count[0:500]:
           if metabolite[1] >= 1000:
               excluded_compounds[metabolite[0]] = metabolite[1]
   
   # print available organisms
   if opts.list_organisms:
       org_list = pythoncyc.all_orgids()
       print clean_ptools_output(org_list)
       exit()
   
   # list of graph edges
   global all_edges
   all_edges = {}
   
   # get all pathways
   pwys = cyc.all_pathways()
   total_pwys = len(pwys)
   pwy_count = 0
   for pwy in pwys:
       print cyc.get_name_string(pwy), pwy_count, "of", total_pwys
       try:
           # get connected components
           connected_components = cyc.pathway_components(pwy)
           for component in connected_components[0]:
               add_edges_from_component(component, pwy)
           exit()
           rxn_list = cyc.get_slot_values(pwy, "REACTION-LIST")
           for i in range(len(rxn_list)):
               for j in range(i+1, len(rxn_list)):
                   # test if pathway reactions should be connected
                   #print rxn_list[i], rxn_list[j], p, "path"
                   if test_for_edge(rxn_list[i], rxn_list[j]):
                       # add edge
                       # print rxn_list[i], rxn_list[j], p, "path"
                       continue
           pwy_count += 1
       except Exception,e:
           print "Warning: Could not get reactions from pathway " + cyc.get_name_string(pwy)
           print str(e)
           exit()
           # rxns = cyc.get_slot_values(p, "REACTION-LIST")
           # for rxn in rxns:
           #     print "LEFT", cyc.get_slot_values(rxn, "LEFT")
           #     print "RIGHT", cyc.get_slot_values(rxn, "RIGHT")
           # continue
   metabolite_count_out = open("/tmp/metabolite_count.pk", "w")
   pickle.dump(metabolite_count, metabolite_count_out)
   metabolite_count_out.close()
   sorted_metabolite_count = sorted(metabolite_count.items(), key=operator.itemgetter(1), reverse=True)
   
   for metabolite in sorted_metabolite_count[0:20]:
       print metabolite
   exit()
   for r in reactions.instances:
       print r.frameid
       left = set(r.left)
       right = set(r.right)
       
       print left.intersection(right)
       


if __name__ == "__main__":
   main()