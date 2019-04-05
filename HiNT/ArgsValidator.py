import os,sys
from HiNT.corelib import *

def opt_validate_hintpre(optparser):
    """Validate arguments from a OptParser object.
    """

    options = optparser.parse_args()
    if options.inputformat == "fastq" and (not options.bwapath or not options.bwaIndex):
        optparser.print_help()
        Info('ERROR: HiNT need the path to bwa and bwa Index to do the alignmet.')
        sys.exit(1)

    if options.inputformat == "cooler" and not options.coolerpath:
        optparser.print_help()
        Info('ERROR: HiNT need the path to cooler to create HI-C contact matrix.')
        sys.exit(1)

    inputdata = options.hicdata
    if options.inputformat == "fastq":
        inputdatas = inputdata.split(',')
        if len(inputdatas) != 2:
            Info('ERROR: Both mates of Hi-C data should be input, a comma ',' should be used to seperate two mates. E.g. hic_1.fastq.gz,hic_2.fastq.gz')
            sys.exit(1)
        else:
            if not os.path.isfile(inputdatas[0]):
                Info('ERROR: %s is not a valid path'%inputdatas[0])
                sys.exit(1)
            if not os.path.isfile(inputdatas[1]):
                Info('ERROR: %s is not a valid path'%inputdatas[1])
                sys.exit(1)
        options.hicdata = inputdatas

    elif options.inputformat == "bam":
        if not os.path.isfile(inputdata):
            Info('ERROR: %s is not a valid path'%inputdata)
            sys.exit(1)
    else:
        Info('ERROR: choose a valid format for input data, fastq or bam')
        sys.exit(1)

    if not options.name:
        options.name = "NA"
        Info("HiNT will use 'NA' as the prefix for all the ouput files")

    if not options.outdir:
        options.outdir = os.getcwd()
        options.outdir = os.path.join(options.outdir,'HiNTpre_OUTPUT')
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)
    else:
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)


    Info("Argument List: ")
    Info("Hi-C data = " + ', '.join(options.hicdata))
    Info("Input format = " + options.inputformat)
    Info("Output format = " + options.outputformat)
    Info("Genome = %s" %options.genome)
    Info("Name = " + options.name)
    Info("Output directory = " + options.outdir)

    return options

def opt_validate_hintcnv(optparser):
    """
    Validate arguments from a OptParser object.
    """

    options = optparser.parse_args()
    if not os.path.isfile(options.matrixfile):
        optparser.print_help()
        Info('ERROR: path to the matrixfile is not valid.')
        sys.exit(1)

    if not os.path.isdir(options.bicseq):
        optparser.print_help()
        Info('ERROR: path to the BICseq package')
        sys.exit(1)

    if not options.name:
        options.name = "NA"
        Info("HiNT will use 'NA' as the prefix for all the ouput files")

    if not options.outdir:
        options.outdir = os.getcwd()
        options.outdir = os.path.join(options.outdir,'HiNTcnv_OUTPUT')
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)
    else:
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)


    Info("Argument List: ")
    Info("Hi-C contact matrix = " + options.matrixfile)
    Info("Hi-C contact matrix format = " + options.format)
    Info("resolution = " + str(options.resolution) + ' kb')
    Info("Genome = " + options.genome)
    Info("BICseq directory = " + options.bicseq)
    Info("Name = " + options.name)
    Info("Output directory = " + options.outdir)

    return options

def opt_validate_hinttransl(optparser):
    """
    Validate arguments from a OptParser object.
    """

    options = optparser.parse_args()
    if options.format == "juicer" and not os.path.isfile(options.matrixfile):
        optparser.print_help()
        Info('ERROR: path to the matrixfile is not valid.')
        sys.exit(1)
    if options.format == "cooler":
        options.matrixfile = options.matrixfile.split(',')
        if len(options.matrixfile) != 2 or not os.path.isfile(options.matrixfile[0]) or not os.path.isfile(options.matrixfile[1]):
            optparser.print_help()
            Info('ERROR: path to the matrixfile is not valid.')
            sys.exit(1)
    if (options.format == "sparse" or options.format == "dense") and not os.path.isdir(options.matrixfile):
        optparser.print_help()
        Info('ERROR: path to the matrixfile is not valid, should be a directory.')
        sys.exit(1)

    if options.chimeric and not os.path.isfile(options.chimeric):
        Info('ERROR: %s is not a valid path'%options.chimeric)
        sys.exit(1)

    if not options.name:
        options.name = "NA"
        Info("HiNT will use 'NA' as the prefix for all the ouput files")

    if not options.outdir:
        options.outdir = os.getcwd()
        options.outdir = os.path.join(options.outdir,'HiNTtransl_OUTPUT')
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)
    else:
        if not os.path.isdir(options.outdir):
            os.mkdir(options.outdir)


    Info("Argument List: ")
    if options.format == "cooler":
        Info("Hi-C contact matrix = " + ', '.join(options.matrixfile))
    if options.format == "juicer":
        Info("Hi-C contact matrix = " + options.matrixfile)
    Info("Hi-C contact matrix format = " + options.format)
    if options.chimeric:
        Info("Hi-C chimeric read pairs = " + options.chimeric)
    Info("Genome = %s" %options.genome)
    Info("Name = " + options.name)
    Info("Output directory = " + options.outdir)

    return options
