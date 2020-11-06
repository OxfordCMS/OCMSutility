'''
combine_lanes.py
=======================

:Author: Nick Ilott
:Tags: Python

Purpose
-------

The purpose of this script is to take a file that contains a mapping between IDs that are given to us from the OGC that relate to fastq filenames and our OCMS IDs and to combine fastq files across multiple lanes such that the resulting fastq files have OCMS identifiers. If the samples have not been run across multiple lanes i.e. there is only a single lane identifier then the script will link files from the original directory to the working directory but with the new identifiers.


Usage
-----

Inputs to the script are a directory containing fastq files that have the WT OGC naming convention and a mapping file (tab delimited) that composes of two columns

+---------+---------+
| WT ID   | OCMS ID |
+---------+---------+

If multiple lanes of sequencing were performed per sample then a new file will be created after concatenating fastq files across lanes. If only one sequencing lane is present per sample then links to the data will be created in the working directory with a new file name as per the mapping file.

Example::

   ocms combine_lanes --fastq-directory=<FASTQ-DIRECTORY> --id-map=<ID_MAPPING_FILE>

Type::

   ocms combine_lanes --help

for command line help.

SDRF output
------------

It is desirable to have metadata information for a sequencing project - this is important as at the end of a project the data are uploaded to a public repository. While all of the metadata cannot be included using this script, it does provide a template if --srdf-output is set to true. This will include minimal information that is required by the EBI ENA for fastq file uploads as well as the files that are available for each sample and the checksums for each file.


Command line options
--------------------

'''

import sys
import os
import glob
import itertools
import cgatcore.experiment as E
import hashlib
import itertools

def readWt2id(wt2id):
    '''
    read identifier map into dictionary
    '''
    mapping = {}
    inf = open(wt2id, encoding="utf_8_sig")
    for line in inf.readlines():
        data = line[:-1].split("\t")
        mapping[data[0]] = data[1]
    return(mapping)

def combineLanes(directory, mapping_dict, read=1):
    '''
    read combination
    '''
    assert read in [1,2], "must specify read pair i.e. 1 or 2"
    if read == 1:
        suffix = "_1.fastq.gz"
        out_suffix = ".fastq.1.gz"
    elif read == 2:
        suffix = "_2.fastq.gz"
        out_suffix = ".fastq.2.gz"

    readfiles = glob.glob(os.path.join(directory, "*"+suffix))
    readfiles.sort()

    found_index = set()
    out_map = open("read"+str(read)+".map", "w")
    for readfile in readfiles:
        prefix = os.path.basename(readfile).replace(suffix, "")
        
        # make sure that the file has a mapping id
        assert prefix in mapping_dict, f"ERROR: File {readfile} exists but There is no sample ID for {prefix} in --id-map file" 
        
        splitname = os.path.basename(readfile).split("_")
        lane, index = splitname[1], splitname[2]
        if index in found_index:
            continue
        else:
            found_index.add(index)
            reads = [x for x in readfiles if index in x]
            tocombine = " ".join(reads)
            newname = mapping_dict[prefix] + out_suffix
            out_map.write(tocombine + " " + newname + "\n")
            if len(reads) == 1:
                assert not os.path.exists(newname), "Remove previously created files from directory. Will not force overwriting"
                statement = f"""ln -s {tocombine}  {newname}"""
            else:
                assert not os.path.exists(newname), "Remove previously created files from directory. Will not force overwriting"
                statement = f"""cat {tocombine} > {newname}"""
            os.system(statement)


def buildSdrf(read1map):
    '''build sdrf template with md5sums for the original
    input files and how they map to sample names
    '''
    infiles = open(read1map)
    
    # peek the number of original files
    noriginal = len(infiles.readline()[:-1].split(" ")) - 1
    noriginal = noriginal*2

    infiles.close()

    # re-open
    infiles = open(read1map)
    
    outf = open("sdrf_template.tsv", "w")
    # write header
    filenumbers = ["File" + str(x) + "\tmd5sums_File" + str(x) for x in range(1, noriginal+1)]
    
    outf.write("\t".join(["Sample",
                         "Organism",
                         "Tissue",
                         "Cell type",
                         "Sex",
                         "Age",
                         "Disease",
                         "Material type",
                         "Library source",
                         "Library selection",
                         "Library layout",
                         "Library orientation",
                         "Library nominal fragment length",
                         "Library nominal fragment Stdev"] + filenumbers) + "\n")

    for line in infiles.readlines():
        data = line[:-1].split(" ")

        # number of original files = no. lanes/runs
        noriginal = len(data)-1
        sample_name = data[-1].replace(".fastq.1.gz", "")

        # files
        orig_files = data[:-1]

        # add second read in pair
        orig_files = orig_files + [x.replace("_1.fastq.gz", "_2.fastq.gz") for x in orig_files]
        hashes = []
        for f in orig_files:
            md5 = hashlib.md5(open(f,'rb').read()).hexdigest()
            hashes.append([os.path.basename(f), md5])
        outf.write("\t".join([sample] + ["NA"]*14 + list(itertools.chain(*hashes))) + "\n")
    outf.close()
            
def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-d", "--fastq-directory", dest="fastq_directory", type=str,
                        help="supply directory containing fastq files to combine")
    parser.add_argument("-i", "--id-map", dest="id_map", type=str,
                        help="supply file containing WT ID to Sample ID (OCMS) mapping")
    parser.add_argument("--srdf-output", dest="sdrf_output", action="store_true",
                        help="flag to specify whether to output a template sdrf file")
    parser.add_argument("--debug", dest="debug", action="store_true",
                        help="turns on full traceback for errors")

    
    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    if not args.debug:
        sys.tracebacklimit=0
            
    wt2id = readWt2id(args.id_map)
    combineLanes(args.fastq_directory, wt2id, read=1)
    combineLanes(args.fastq_directory, wt2id, read=2)

    if args.sdrf_output:
        E.warn("Building sdrf_template.tsv: This may take some time :)")
        buildSdrf("read1.map")
    
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))



