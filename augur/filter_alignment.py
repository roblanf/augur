"""
Align multiple sequences from FASTA.

Input: a FASTA alignment
Output: a path to a filtered alignment with sites masked by N's and strains removed

"""

import os, glob, Bio
from shutil import copyfile, rmtree
import numpy as np
from Bio import AlignIO, SeqIO, Seq, Align
from .utils import run_shell_command, nthreads_value, shquote
from .tree import load_excluded_sites, mask_sites_in_multiple_sequence_alignment
from .align import read_alignment
from . import logtools

log = logtools.get_logger()


class AlignmentError(Exception):
    # TODO: this exception should potentially be renamed and made augur-wide
    # thus allowing any module to raise it and have the message printed & augur
    # exit with code 1
    pass

def register_arguments(parser):
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--alpha', type=float, default=0.05, help="(default is 0.05), exclude sequences that significantly increase the tree diameter, with a false-positive rate of alpha. This parameter is equivalent to the -q parameter in TreeShrink. Set this to <0 (e.g. -1) if you don't want to run TreeShrink")
    parser.add_argument('--gapthresh', type=str, default=0.2, help='only keep columns with <= this fraction of gaps [default 0.2], using EASEL alimask. Set to 1.0 to mask no columns with gaps.')
    parser.add_argument('--minlength', type=str, default=1, help='remove strains with a sequence length < this, using EASEL alimanip. Leave at the default value of 1 to remove no strains based on length.')
    parser.add_argument('--maxambig', type=str, default=None, help='remove strains with >= this many ambiguous residues, using EASEL alimanip. The default (None) means no strains will be filtered based on ambiguities')
    parser.add_argument('--exclude_strains', type=str, help="file with list of strains that are to be excluded regardless of other filtering steps. One strain per line.")
    parser.add_argument('--retain_strains', type=str, help="file with list of strains that are to be retained regardless of other filtering steps. One strain per line.")
    parser.add_argument('--exclude_sites', type=str, help="file with list of sites that are to be excluded regardless of other filtering steps. This should be one of: a BED file (with a .bed extension), a tab-delimited DRM file, or a plain text file with one position per line (1-indexed in the last case).")
    parser.add_argument('--retain_sites', type=str, help="file with list of sites (1-referenced relative to the input alignment) that are to be retained regardless of other filtering steps. This should be one of: a BED file (with a .bed extension), a tab-delimited DRM file, or a plain text file with one position per line (1-indexed in the last case).")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--output', '-o', help="output file", required=True)


def run(args):
    '''
    filter an alignment
    '''

    log.info("Starting alignment filtering")


    # Convert  VCF to fasta, for now this just uses identical code to that in tree.py
    # TODO: break out the below section as a function and put it in utils.py
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):

        log.info("VCF alignment detected")

        # Prepare a multiple sequence alignment from the given variants VCF and
        # reference FASTA.
        if not args.vcf_reference:
            log.error("a reference Fasta is required with VCF-format alignments, please supply one with --vcf-reference")
            return 1
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
        aln = sequences
    else:
        # Use the multiple sequence alignment as is.
        aln = args.alignment

    # construct reduced alignment if needed
    if is_vcf:
        log.info("Converting VCF alignment to FASTA of variable sites")
        variable_fasta = write_out_informative_fasta(compress_seq, args.alignment, stripFile=args.exclude_sites)
        fasta = variable_fasta
    else:
        fasta = aln

    #############
    # FILTERING #
    #############

    all_aln = read_alignment(fasta)

    if args.maxambig == None:
        args.maxambig = all_aln.get_alignment_length()

    # run alimanip to remove strains based on cutoffs supplied by user
    with(logtools.indented(log, "Removing strains based on minimum length %s and maximum number of ambiguities %s" %(args.minlength, args.maxambig))):

        fasta_alimanip = "".join([fasta, "_alimanip.fasta"])

        call = ["esl-alimanip", "--lmin", str(args.minlength), 
                                "--xambig", str(args.maxambig),
                                "--informat afa", 
                                "--outformat afa", 
                                "--dna", 
                                "-o", fasta_alimanip, 
                                fasta]

        cmd = " ".join(call)
        log.info("Running esl-alimask with the following command:")
        log.info(cmd)
        log.info("citation for esl-alimask: https://github.com/EddyRivasLab/easel")

        run_shell_command(cmd, raise_errors = True)

        # make an alignment that removes taxa suggested by alimanip and user excluded strains, but retains user retained clades
        all_names = [x.description for x in all_aln]
        ali_names = [x.description for x in read_alignment(fasta_alimanip)]
        ali_remove = set(all_names) - set(ali_names)

        log.info("alimanip indentified %d strains to exclude based on minimum length and/or maximum number of ambiguities" % len(ali_remove))

        # don't remove those the user wants to keep
        if args.retain_strains:
            log.debug("loading user_retain_strains file")
            with open(args.retain_strains) as f:
                user_retain_strains = [x.strip() for x in f]
        else:
            user_retain_strains = []

        if args.exclude_strains:
            log.debug("loading user_exclude_strains file")
            with open(args.exclude_strains) as f:
                user_exclude_strains = [x.strip() for x in f] 
        else:
            user_exclude_strains = []

        remove = ali_remove - set(user_retain_strains)

        log.info("user supplied %d strains to exclude and %s strains to retain", 
            len(user_exclude_strains), 
            len(user_retain_strains))

        
        log.info("combining all information, %d (%.2f percent) strains of the total %d will be removed",
            len(remove),
            100*float(len(remove)) / len(all_names),
            len(all_names))

        # add to the list seqs user wants to exclude (but not if they're not in the alignment to start with)
        for strain in user_exclude_strains:
            if strain in all_names:
                remove.add(strain)

        keep = set(all_names) - remove
        
        # write out a new alignment with these strains excluded
        log.debug("writing alimanip filtered alignment in FASTA to %s", fasta_alimanip)

        with open(fasta_alimanip, "w") as oh:
            for record in all_aln:
                if record.description in keep:
                    Bio.SeqIO.write(record, oh, "fasta")

        removed_strains_file = "".join([fasta, "-excluded_strains.txt"])
        log.debug("writing excluded strains file to %s", removed_strains_file)
        with open(removed_strains_file, 'w') as f:
            f.write("strain\treason\n")
            for s in remove:
                if s in user_remove_strains:
                    f.write("\t".join([s, "excluded by user"]))
                else:
                    f.write("\t".join([s, "'length < %s and/or number of ambiguities > %s filtered with esl-alimanip'\n" %(str(args.minlength), str(args.maxambig))]))

        log.info("file with information on each strain excluded is here: %s", removed_strains_file)

    # Run EASEL, make a set of sites to be deleted from input alignment. 1 for exclude, 0 for include
    with(logtools.indented(log, "Removing sites with a proportion of gaps greater than %s" % str(args.gapthresh))):

        # note we pipe to nowhere because esl-alimask streams the output alignment, which we don't want in this case
        esl_maskfile = "".join([fasta_alimanip, "-eslMask"])
        call = ["esl-alimask", "--gapthresh", str(args.gapthresh), "--informat afa", "--dna", "--fmask-all", esl_maskfile, "-g", fasta_alimanip, "> /dev/null 2>&1"]
        cmd = " ".join(call)
        log.info("Running esl-alimanip with the following command:")
        log.info(cmd)
        log.info("citation for esl-alimanip: https://github.com/EddyRivasLab/easel")

        run_shell_command(cmd, raise_errors = True)
        log.debug("loading esl-alimanip mask file from %s", esl_maskfile)
        with open(esl_maskfile) as f:
            esl_mask = f.read()
        # clean up and make the 1's mean 'remove' and the zero's mean 'keep'
        esl_mask = esl_mask.rstrip("\n")
        esl_mask = esl_mask.replace("0", "f")
        esl_mask = esl_mask.replace("1", "0")
        esl_mask = esl_mask.replace("f", "1")
        esl_maskarray = np.array(list(map(int, list(esl_mask))))
        os.remove(esl_maskfile) # clean up

        maskarray = esl_maskarray
        sites_to_remove = np.unique(np.where(maskarray>0))

        log.info("esl-alimanip found %d sites with a proportion of gaps greater than %s",
            len(sites_to_remove),
            str(args.gapthresh))


        log.debug("loading user-supplied exclude_sites file from %s" %(args.exclude_sites))
        # add sites the user wants to exclude 
        if args.exclude_sites:
            user_exclude_sites = load_excluded_sites(args.exclude_sites)
            sites_to_remove = np.unique(np.append(sites_to_remove, user_exclude_sites))
        else:
            user_exclude_sites = []
        
        log.debug("loading user-supplied retain_sites file from %s" %(args.retain_sites))
        # remove from the list the sites the user wants to include    
        if args.retain_sites:
            user_retain_sites = load_excluded_sites(args.retain_sites)
            sites_to_remove = np.delete(sites_to_remove, np.where(np.in1d(sites_to_remove, user_retain_sites)))
        else:
            user_retain_sites = []

        log.info("user supplied %d sites to exclude and %s sites to retain", 
            len(user_exclude_sites), 
            len(user_retain_sites))

        log.info("combining all information, %d (%.2f percent) sites of the total %d will be masked",
            len(sites_to_remove),
            100*float(len(sites_to_remove)) / len(esl_maskarray),
            len(esl_maskarray))

        log.debug("masking sites from alignment with biopython")

        if len(sites_to_remove>0):
            masked_sites_file_user = "".join([fasta, "-masked_sites.txt"])
            
            log.info("file with a list of all sites masked is here: %s", masked_sites_file_user)
            
            sites_to_remove_user = sites_to_remove + 1
            np.savetxt(masked_sites_file_user, sites_to_remove_user.astype(int), fmt='%i')

            # we also need to save a list of 0 indexed sites for masking the alignment in augur
            masked_sites_file = "".join([fasta, "-masked_sites_0based.txt"])
            np.savetxt(masked_sites_file, sites_to_remove.astype(int), fmt='%i')

            # this alignment is masked with N's
            masked_sites_aln = mask_sites_in_multiple_sequence_alignment(fasta_alimanip, masked_sites_file)

        else:
            masked_sites_aln = fasta_alimanip # nothing to do here


    # Run TreeShrink on the masked alignment, make a list of taxa to be deleted from input alignment
    if(args.alpha>0):
        with(logtools.indented(log, "Removing taxa from masked alignment using TreeShrink to detect weirdly long branches (alpha=%s)" % str(args.alpha))):

            # make a quick tree with IQ-TREE, designed to operate very much like fastME here
            log.info("building fast BIONJ tree in IQ-TREE using HKY model: good enough to detect long branches")
            call = ["iqtree", "-s", masked_sites_aln, "-te BIONJ", "-m HKY", "-nt", str(args.nthreads), "-redo"]
            cmd = " ".join(call)
            run_shell_command(cmd, raise_errors = True)
            bionj_tree = "".join([masked_sites_aln, ".treefile"])

            ts_outdir =  os.path.join(os.path.dirname(bionj_tree), "".join([os.path.basename(bionj_tree), "_treeshrink"]))
            call = ["run_treeshrink.py", "-t", bionj_tree, "-q", str(args.alpha), "-c", "-o", ts_outdir]
            cmd = " ".join(call)

            log.info("Running TreeShrink with the following command:")
            log.info(cmd)
            log.info("citation for TreeShrink: Mai, Uyen, and Siavash Mirarab. TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees. BMC genomics 19.5 (2018): 272.")

            run_shell_command(cmd, raise_errors = True)

            ts_maskfile = glob.glob(os.path.join(ts_outdir, "*.txt"))[0]
            with open(ts_maskfile) as f:
                ts_mask = f.read()
            ts_mask = ts_mask.rstrip('\n')
            ts_mask = ts_mask.rstrip('\t')

            ts_remove = set(ts_mask.split("\t"))

            # Never exclude taxa the user wants to keep
            remove = ts_remove - set(user_retain_strains)

            all_aln = read_alignment(masked_sites_aln)

            all_names = [x.description for x in all_aln]

            keep = set(all_names) - remove

            # write out a new alignment with these strains excluded
            with open(args.output, "w") as oh:
                for record in all_aln:
                    if record.description in keep:
                        Bio.SeqIO.write(record, oh, "fasta")

            log.info("TreeShrink indentified %d strains to exclude based on inducing weirdly long branches" % len(ts_remove))
            log.info("Combined with information on strains user wants to retain, %d strains will be removed" % len(remove))


            with open(removed_strains_file, 'a') as f:
                for s in remove:
                    f.write("\t".join([s, "'filtering with TreeShrink'\n"]))
            log.info("information on strains removed by TreeShrink has been added to: %s", removed_strains_file)

    else:
        log.info("Not running TreeShrink because alpha was set to %s" % str(args.alpha))
        os.rename(masked_sites_aln, args.output)

    # clean up
    log.debug("cleaning up by deleting pointless files")
    os.remove(fasta_alimanip)
    os.remove(masked_sites_file)
    junk = glob.glob(os.path.join(os.path.dirname(fasta), "masked_*"))
    for name in junk:
        if name != args.output:
            if name.endswith("alimanip.fasta.treefile_treeshrink"):
                rmtree(name)
            else:
                os.remove(name)

    log.info("Alignment filtering complete")
    log.info("Alignment with strains removed and sites masked is here: %s" %(args.output))
    return args.output 