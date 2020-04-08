"""
Align multiple sequences from FASTA.

Input: a FASTA alignment
Output: a path to a filtered alignment with sites masked by N's and strains removed

"""

import os, glob, Bio
from shutil import copyfile
import numpy as np
from Bio import AlignIO, SeqIO, Seq, Align
from .utils import run_shell_command, nthreads_value, shquote
from .tree import load_excluded_sites
from .align import read_alignment

class AlignmentError(Exception):
    # TODO: this exception should potentially be renamed and made augur-wide
    # thus allowing any module to raise it and have the message printed & augur
    # exit with code 1
    pass

def register_arguments(parser):
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--filters', default='all', choices=["treeshrink", "gblocks", "easel", "all"], help="tree builder to use")
    parser.add_argument('--alpha', type=float, default=0.05, help="(default is 0.05), exclude sequences that significantly increase the tree diameter, with a false-positive rate of alpha. This will be done on a BIONJ tree created with IQ-TREE, or on the tree supplied with --tree. This parameter is equivalent to the -q parameter in TreeShrink")
    parser.add_argument('--gapthresh', type=str, default=0.2, help='only keep columns with <= <x> fraction of gaps in them [default 0.2], using EASEL alimask.')
    parser.add_argument('--lmin', type=str, default=1, help='remove strains w/length < <n> residues using EASEL alimanip (done after all other operations).')
    parser.add_argument('--xambig', type=str, default=100000000, help='remove strains with >= <n> ambiguous residues using EASEL alimanip (done after all other operations).')
    parser.add_argument('--treat_n_as_gap', type=str, default=True, help="Treat N's as gaps? Recommended since ML tree inference treats both the same")
    parser.add_argument('--exclude_strains', type=str, help="file with list of strains that are to be excluded regardless of other filtering steps. One strain per line.")
    parser.add_argument('--retain_strains', type=str, help="file with list of strains that are to be retained regardless of other filtering steps. One strain per line.")
    parser.add_argument('--exclude_sites', type=str, help="file with list of sites that are to be excluded regardless of other filtering steps. This should be one of: a BED file (with a .bed extension), a tab-delimited DRM file, or a plain text file with one position per line (1-indexed in the last case).")
    parser.add_argument('--retain_sites', type=str, help="file with list of sites (1-referenced relative to the input alignment) that are to be retained regardless of other filtering steps. This should be one of: a BED file (with a .bed extension), a tab-delimited DRM file, or a plain text file with one position per line (1-indexed in the last case).")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                                help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--gblocks-args', type=str, default="", help='extra arguments to be passed directly to the gblocks, which otherwise runs with default parameters except for allowing all gaps (since this is treated later with the --gapthresh parameter) (e.g., --gblocks-args="-b4=20")')
    parser.add_argument('--output', '-o', help="output file", required=True)


def run(args):
    '''
    filter an alignment
    '''

    # Convert  VCF to fasta, for now this just uses identical code to that in tree.py
    # TODO: break out the below section as a function and put it in utils.py
    # check alignment type, set flags, read in if VCF
    is_vcf = False
    ref = None
    if any([args.alignment.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        # Prepare a multiple sequence alignment from the given variants VCF and
        # reference FASTA.
        if not args.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return 1
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        is_vcf = True
        aln = sequences
    elif args.exclude_sites:
        # Mask excluded sites from the given multiple sequence alignment.
        aln = mask_sites_in_multiple_sequence_alignment(args.alignment, args.exclude_sites)
    else:
        # Use the multiple sequence alignment as is.
        aln = args.alignment

    if args.output:
        tree_fname = args.output
    else:
        tree_fname = '.'.join(args.alignment.split('.')[:-1]) + '.nwk'

    # construct reduced alignment if needed
    if is_vcf:
        variable_fasta = write_out_informative_fasta(compress_seq, args.alignment, stripFile=args.exclude_sites)
        fasta = variable_fasta
    else:
        fasta = aln

    #############
    # FILTERING #
    #############

    # run alimanip to remove strains based on cutoffs supplied by user
    if args.filters in ['all', 'easel']:

        fasta_alimanip = "".join([fasta, "_alimanip.fasta"])

        call = ["esl-alimanip", "--lmin", str(args.lmin), 
                                "--xambig", str(args.xambig),
                                "--informat afa", 
                                "--outformat afa", 
                                "--dna", 
                                "-o", fasta_alimanip, 
                                fasta]

        cmd = " ".join(call)
        run_shell_command(cmd, raise_errors = True)

        # make an alignment that removes taxa suggested by alimanip and user excluded strains, but retains user retained clades
        all_aln = read_alignment(fasta)
        all_names = [x.description for x in all_aln]
        ali_names = [x.description for x in read_alignment(fasta_alimanip)]
        ali_remove = set(all_names) - set(ali_names)

        # don't remove those the user wants to keep
        if args.retain_strains:
            with open(args.retain_strains) as f:
                user_retain_strains = [x.strip() for x in f]
        else:
            user_retain_strains = []

        if args.exclude_strains:
            with open(args.exclude_strains) as f:
                user_exclude_strains = [x.strip() for x in f] 
        else:
            user_exclude_strains = []

        remove = ali_remove - set(user_retain_strains)

        # add to the list seqs user wants to exclude (but not if they're not in the alignment to start with)
        for strain in user_exclude_strains:
            if strain in all_names:
                remove.add(strain)

        keep = set(all_names) - remove

        # write out a new alignment with these strains excluded
        with open(fasta_alimanip, "w") as oh:
            for record in all_aln:
                if record.description in keep:
                    Bio.SeqIO.write(record, oh, "fasta")

        removed_strains_file = "".join([fasta, "_removed_strains.txt"])
        with open(removed_strains_file, 'w') as f:
            f.write("strain\treason\n")
            for s in remove:
                f.write("\t".join([s, "'lmin and xambig filtering with esl-alimanip'\n"]))


        print("%d strains removed using lmin and xambig filtering: " % len(remove))
        print("file with information on strains removed is here: ", removed_strains_file)

    # Run gblocks, make a set of sites to be deleted from input alignment
    if args.filters in ['all', 'gblocks']:
        call = ["gblocks", fasta_alimanip, "-k=y", "-t=d", "-p=n", "-s=y", "-b5=a", args.gblocks_args]
        cmd = " ".join(call)
        print("Running gblocks via:\n\t" + cmd +
              "\n\tCastresana, J. (2000). Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis." +
              "\n\tMolecular Biology and Evolution 17, 540-552.\n")

        # set raise_errors to False because gblocks seems to trigger run_shell_command to report an error when it runs fine. 
        run_shell_command(cmd, raise_errors = False)


        # extract sequence P1;Gblocks, which is a mask with .==remove and #==retain
        # convert to .=1 and #=0 and make it a vector of ints
        gb_maskfile = "".join([fasta_alimanip, "-gbMask"])
        seqs = read_alignment(gb_maskfile)
        gb_mask = str(seqs[-1].seq)[:-1] # mask is the last seq in the fasta, and gblocks adds a "*" that we have to remove
        
        print(gb_mask)
        print(len(gb_mask))
        print(seqs)
        egg = read_alignment(fasta_alimanip)
        print(egg)

        gb_mask = gb_mask.replace("#", "0")
        gb_mask = gb_mask.replace(".", "1")
        gb_maskarray = np.array(list(map(int, list(gb_mask))))
        #os.remove(gb_maskfile) # clean up
    else:
        gb_maskarray = np.empty(1)

    # Run EASEL, make another set of sites to be deleted from input alignment. as above, 1 for exclude, 0 for include
    if args.filters in ['all', 'easel']:
        # note we pipe to nowhere because esl-alimask streams the output alignment, which we don't want in this case
        esl_maskfile = "".join([fasta_alimanip, "-eslMask"])
        call = ["esl-alimask", "--gapthresh", str(args.gapthresh), "--informat afa", "--dna", "--fmask-all", esl_maskfile, "-g", fasta_alimanip, "> /dev/null 2>&1"]
        cmd = " ".join(call)
        print("Running esl-alimask via:\n\t" + cmd +
              "\n\thttps://github.com/EddyRivasLab/easel\n")
        run_shell_command(cmd, raise_errors = True)
        with open(esl_maskfile) as f:
            esl_mask = f.read()
        # clean up and make the 1's mean 'remove' and the zero's mean 'keep'
        esl_mask = esl_mask.rstrip("\n")
        esl_mask = esl_mask.replace("0", "f")
        esl_mask = esl_mask.replace("1", "0")
        esl_mask = esl_mask.replace("f", "1")
        esl_maskarray = np.array(list(map(int, list(esl_mask))))
        #os.remove(esl_maskfile) # clean up

    else:
        esl_maskarray = np.empty(1)

    # sum the two exclude lists, such that any value >0 means some method thought we should exclude it
    maskarray = gb_maskarray + esl_maskarray
    sites_to_remove = np.unique(np.where(maskarray>0))

    # add sites the user wants to exclude 
    if args.exclude_sites:
        excl_sites = load_excluded_sites(args.exclude_sites)
        sites_to_remove = np.unique(np.append(sites_to_remove, excl_sites))
    
    # remove from the list the sites the user wants to include    
    if args.include_sites:
        incl_sites = load_excluded_sites(args.include_sites)
        sites_to_remove = np.delete(sites_to_remove, np.where(np.in1d(sites_to_remove, incl_sites)))

    # we then make the cut alignment
    print("removing %d sites (%.2f percent) from the input alignment" % (len(sites_to_remove), (float(len(sites_to_remove) / len(esl_maskarray)))))

    if len(sites_to_remove>0):
        masked_sites_file_user = "".join([fasta, "-filtered_sites.txt"])
        print("list of sites removed by filtering is here:", masked_sites_file)
        sites_to_remove_user = sites_to_remove + 1
        np.savetxt(masked_sites_file_user, sites_to_remove_user.astype(int), fmt='%i')

        # we also need to save a list of 0 indexed sites for masking the alignment in augur
        masked_sites_file = "".join([fasta, "-filtered_sites_0based.txt"])
        np.savetxt(masked_sites_file, sites_to_remove.astype(int), fmt='%i')

        # this alignment is masked with N's
        masked_sites_aln = mask_sites_in_multiple_sequence_alignment(fasta_alimanip, masked_sites_file)

    else:
        masked_sites_aln = fasta_alimanip # nothing to do here






    # Run TreeShrink on the masked alignment, make a list of taxa to be deleted from input alignment
    if args.filters in ['all', 'treeshrink']:

        # make a quick tree with IQ-TREE, designed to operate very much like fastME here
        call = ["iqtree", "-s", masked_sites_aln, "-te BIONJ", "-m HKY", "-nt", args.nthreads]
        cmd = " ".join(call)
        run_shell_command(cmd, raise_errors = True)
        bionj_tree = "".join([masked_sites_aln, ".treefile"])

        print("Running TreeShrink via:\n\t" + cmd +
            "\n\tMai, Uyen, and Siavash Mirarab. TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees." +
            "\n\tBMC genomics 19.5 (2018): 272.\n")
        ts_outdir =  os.path.join(os.path.dirname(bionj_tree), "".join([os.path.basename(bionj_tree), "_treeshrink"]))
        call = ["run_treeshrink.py", "-t", bionj_tree, "-q", str(alpha), "-c", "-o", ts_outdir]
        cmd = " ".join(call)
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

        print("%d strains removed using TreeShrink filtering: " % len(remove))
        print("file with information on strains removed is here: ", removed_strains_file)

        with open(removed_strains_file, 'a') as f:
            f.write("strain\treason\n")
            for s in remove:
                f.write("\t".join([s, "'lmin and xambig filtering with esl-alimanip'\n"]))

    else:
        
        os.rename(masked_sites_aln, args.output)


    return final_align 