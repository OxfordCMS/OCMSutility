'''
ocms_diff2graphlan.py
=============================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The purpose of this script is to take results of a differential abundance test and Create graphical parameters for input to graphlan. This allows for differential features to be visualised down the taxonomic tree.


Inputs
-------

The input to ocms_diff2graphlan.py is a tab-delimited text file of differential abundance results. These results should have been generated for each taxonomic level of interest, concantenated and have the following column headers:


+-----------+-----------+----------------+-----------+----------+----------+
| taxa      | baseMean  | log2FoldChange | lfcSE     | pvalue   | padj     |
+-----------+-----------+----------------+-----------+----------+----------+

The reason for these headers is that the script was based on the output from DESeq2. If you have used something else for differential abundance testing then you can have columns with NA values for most of these columns. The only ones used by the script are log2FoldChange and padj. Names must be the same though.

The taxa column has to have taxa names in a specific format which is based on the format that defines a tree structure e.g. Bacteria.Firmicutes.Clostridia represents the class Clostridia. 


Ouputs
-------

The main output is an tab-delimited text file that will serve as input to graphlan. It also outputs an input_tree.txt file which is the other input into graphlan_annotate.py.

Usage
-----

This is an example of how to use ocms_diff2graphlan.py. This example will run through the generation of the annotations and input_tree.txt and how these are then input into graphlan_annotate.py and finally graphlan.py for the visualisation. Things are a little tricky because graphlan is only compatible with python-2.7. This means that after ocms_diff2graphlan.py has been run you will need to switch to an environment that allows you to run graphlan.

You can use the test data in this repository to test this out.


Example::

   /gfs/devel/nilott/OCMS_Sandbox/Py_utility/scripts/ocms_diff2graphlan.py --diff-table=diff_table.tsv --tax-tree=tree.txt --highest-level=phylum -u padj > annot.txt

This is taking the differential abundance table and taxonomic tree inputs and will output annotations for all clades up to phylum level. Colours and bars in teh annotated tree will be based on the adjusted p-value in the diff-table i.e. red for significantly increased and blue for significantly decreased. These annotations are written to annot.txt that will be the input into graphlan_annotate.py.


In order to run graphlan you will need to set up a python-2.7 virtual environment. After deactivating any conda environments that you have loaded do::

    module load apps/all; module load apps/python/2.7.13

Create a virtual environment::

    virtualenv graphlan_analysis

Then you can activate this::

    source graphlan_analysis/bin/activate

and install graphlan::

    pip install graphlan

Once this is installed you can use graphlan to annotate your differential abundance tree::

    graphlan_annotate.py --annot annot.txt input_tree.txt output_tree.xml 

This will produce the .xml file that graphlan will use to draw the tree with annotations::

    graphlan.py output_tree.xml annotated_tree.png


This will give you hopefully the desired output of a pretty annotated tree.





Type::

   python diff2graphlan_annotations.py --help

for command line help.

Command line options
--------------------

'''

import sys
import math
import cgatcore.iotools as IOTools
import collections
import cgatcore.experiment as E
from matplotlib import colors
import six
import random
import string

def readTree(infile, highest_level="family"):
    '''
    read the tree file and map nodes
    to use specified highest node.
    '''
    taxa = collections.defaultdict(set)
    inf = IOTools.open_file(infile)
    inf.readline()
    E.info("Using taxon names from --tax-tree to build taxonomic tree")    
    for line in inf.readlines():

        # the data have to contain at least kingdom
        # check this here
        assert "k__" in line, "malformatted data, data must start at kingom"
        
        data = line.strip("\n")
        
        # use the length of the taxon name
        # to determine the taxonomic level
        # and add additional lower ranks
        # depending if they contain this string
        name_length = len(data.split("."))
        if name_length == 1:
            level = "kingdom"
        elif name_length == 2:
            level = "phylum"
        elif name_length == 3:
            level = "class"
        elif name_length == 4:
            level = "order"
        elif name_length == 5:
            level = "family"
        elif name_length == 6:
            level = "genus"
        elif name_length == 7:
            level = "species"
        else:
            try:
                line = line.replace("._", "_")
            except ValueError("malformatted line %s" % line):
                continue

        level_indices = {"kingdom": 0, "phylum": 1, "class": 2, "order": 3, "family": 4, "genus": 5, "species": 6}
        highest_level_index = level_indices[highest_level]
        if level_indices[level] == highest_level_index:
            taxa[data].add(data)
        elif level_indices[level] > highest_level_index:
            taxon = ".".join(data.split(".")[0:highest_level_index+1])
            taxa[taxon].add(data)
        else:
            continue
    return(taxa)
        

#############################################################
#############################################################
#############################################################

def getColours(ncols):
    '''
    get colours - lower clades inherit from "highest level"
    clade
    '''
    colors_ = list(six.iteritems(colors.cnames))

    # Add the single letter colors.
    for name, rgb in six.iteritems(colors.ColorConverter.colors):
        hex_ = colors.rgb2hex(rgb)
        colors_.append((name, hex_))

    # Transform to hex color values.
    hex_ = [color[1] for color in colors_]
    cols = hex_[0:ncols]

    return cols
    
#############################################################
#############################################################
#############################################################
    

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-u", "--use", dest="use", type=str,
                      choices=("pval", "padj"),
                      help="Type of p-value to use for clade size")

    parser.add_argument("-d", "--diff-table", dest="diff_table", type=str,
                      help="differential abundance table")

    parser.add_argument("-t", "--tax-tree", dest="tax_tree", type=str,
                      help="full taxonomic tree")

    parser.add_argument("-l", "--highest-level", dest="highest_level", type=str,
                      help="highest taxonomic level to visualise")

    parser.add_argument("-f", "--filter", dest="filter", action="store_true",
                      help="do you want to filter? will filter based on highest-level")

    parser.add_argument("-k", "--keep", dest="keep", type=str,
                      help="keep all clades below these")

    parser.add_argument("--additional-labels", dest="additional_labels", type=str,
                      help="by default just the highest level labels are shown. Here you can add additional labels")
    
    
    # add common args (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    # read tree
    tree = readTree(args.tax_tree, args.highest_level)

    # filter if neccessary
    new_tree = {}
    keep = set()
    if args.filter:
        assert args.keep, "must specify which clades to keep"
        to_keep = args.keep.split(",")
        for t in to_keep:
            new_tree[t] = tree[t]
    else:
        assert len(list(tree.keys())) < 159, "not enough colours to support n = %i clades, please filter" % len(tree)
        new_tree = tree

    #get those to keep
    keep = set()
    for taxon, rest in new_tree.items():
        keep.add(taxon)
        for r in rest:
            keep.add(r)

    # get colours
    ncols = len(list(new_tree.keys()))
    colours = getColours(ncols)
    taxon2colour = {}
    for i in range(ncols):
        taxon2colour[list(new_tree.keys())[i]] = colours[i]

    # add in colours for all clade nodes that is
    # based on the highest node
    for h, r in new_tree.items():
        for taxon in r:
            taxon2colour[taxon] = taxon2colour[h]

    # read diff and output annotations
    result = collections.defaultdict(list)
    ps = []
    fcs =[]
    taxa = []
    colours = []
    shapes = []
    sig = []
    for line in open(args.diff_table).readlines():

        # skip header
        if "taxa" and "log2FoldChange" in line:
            continue

        data = line.strip("\n").split("\t")
        taxon = data[0]
        if taxon not in keep:
            continue

        # assign colour
        if taxon in list(taxon2colour.keys()):
            colour = taxon2colour[taxon]
        else:
            colour = "NA"
        colours.append(colour)
            
        # append taxon to list
        taxa.append(taxon)
        
        # use -log10 p-value for clade size
        if args.use == "pval":
            column = 5
        elif args.use == "padj":
            column = 6
        else:
            raise ValueError("must use pval or padj, not %s for clade size" % args.use)

        # catch NA pvalues
        if data[column] == "NA":
            data[column] = 1
        
        p = -math.log10(float(data[column]))
        if p >= 1.3:
            shapes.append("*")
            sig.append(taxon)
        else:
            shapes.append("o")

        p = p*100
        result[taxon].append(p)
        ps.append(p)

        # fold changes
        if data[2] == "NA":
            fc = 0
        else:
            fc = data[2]
        fcs.append(fc)
        
    # output annotations
    args.stdout.write("%s\t%s\n" % ("clade_separation", "0.9"))
    for t, s in zip(taxa, shapes):
         args.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_shape", s))
    for t, c in taxon2colour.items():
         args.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", c))
    for t, c in taxon2colour.items():
        args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", c))
    for t, c in taxon2colour.items():
        args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_font_size", 12))

    for t, p, f in zip(taxa, ps, fcs):
        if t in sig and float(f) > 0:
            args.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", "r"))
            args.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", 200))
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", "r"))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_height", 1, f))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_color", 1, "r"))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_alpha", 1, 0.5))

        elif t in sig and float(f) < 0:
            args.stdout.write("%s\t%s\t%s\n" % (t, "clade_marker_color", "b"))
            args.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", 200))
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_background_color", "b"))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_height", 1, float(f)*-1))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_color", 1, "b"))
            args.stdout.write("%s\t%s\t%s\t%s\n" % (t, "ring_alpha", 1, 0.5))
            
        elif t not in sig:
            args.stdout.write("%s\t%s\t%f\n" % (t, "clade_marker_size", p))

    # only output annotation for highest-level and
    # additional labels
    if args.additional_labels:
        additional_labels = args.additional_labels.split(",")
    else:
        additional_labels = []
    for t, p in zip(taxa, ps):
        if t in additional_labels:

            # Plans to make label keys for some things rather than full names
#           if "_" in t:
#               a =  "".join(random.sample(list(string.ascii_lowercase),2)) + ":" + t.split(".")[-1]
#           else:
            a = t.split(".")[-1]
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation", a))
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_rotation", 90))
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_font_size", 8))
        elif t in list(tree.keys()):
            a = t.split(".")[-1]
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation", a))
            args.stdout.write("%s\t%s\t%s\n" % (t, "annotation_font_size", 12))
            
    # write the tree out
    outf = open("input_tree.txt", "w")
    for x, y in new_tree.items():
        for taxon in y:
            outf.write(taxon + "\n")
    outf.close()
           
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
