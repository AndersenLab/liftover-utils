# Liftover 
import sys, os
import tempfile
import subprocess
from subprocess import *

# Check to see if CHROM DIFFs are available.
if os.path.isfile("remap_gff_between_releases.pl") == False:
    os.system("wget ftp://ftp.sanger.ac.uk/pub2/wormbase/software/Remap-between-versions/remap.tar.bz2 && gunzip -f remap.tar.bz2 && tar -xf remap.tar")
    os.system("mv Remap-for-other-groups/remap_gff_between_releases.pl remap_gff_between_releases.pl")
    os.system("mv Remap-for-other-groups/CHROMOSOME_DIFFERENCES/ CHROMOSOME_DIFFERENCES/")
    os.system("rm -f -r Remap-for-other-groups/")
    os.remove("remap.tar")

if len(sys.argv) == 1:
    print """

        usage:

            liftover.py <file> <release1> <release2> <chrom_col> <pos_start_column> [pos_end_column]

            OR

            liftover.py <file> <release1> <release2> {BCF|VCF|GFF} 

        Assumes '+' strandedness

    """
elif len(sys.argv) == 5:
    if (sys.argv[4] in ["BCF", "VCF", "GFF"]):
        vcf = True
        chrom_col, start_col, end_col = 0, 1, 1
        delim = "\t"
        os.system("bcftools query -f '%%CHROM\t%%POS\n' %s > bcf_pos.txt" % sys.argv[1])
        variant_positions = file("bcf_pos.txt",'r')
    else:
        raise Exception("You must specify a valid file type: VCF or GFF.")
        
else:
    pass


def pipe_out(line):
    try:
        sys.stdout.write(line + "\n")
    except IOError:
        try:
            sys.stdout.close()
        except IOError:
            pass
        try:
            sys.stderr.close()
        except IOError:
            pass

gff_temp = tempfile.NamedTemporaryFile().name
gff_liftover = tempfile.NamedTemporaryFile().name
gff = file(gff_temp, 'w+')
release1, release2 = sys.argv[2:4]


for l in variant_positions.xreadlines():
    if l.startswith("#") == False:
        l = l.replace("\n","").split(delim)
        if l[0].lower() == "chrm":
            l[0] = "CHROMOSOME_MtDNA"
        # Write out the coordinates in temporary gff file.
        line_out = "%s\t.\t.\t%s\t%s\t.\t+\t.\t%s\t%s\t%s\n" % tuple([l[chrom_col], l[start_col], l[end_col]]*2)
        gff.write(line_out)

gff.close()


# Generate Liftover Coordinates
remap_command = "perl remap_gff_between_releases.pl -gff=%s -release1=%s -release2=%s -output=%s" % (gff_temp, release1, release2, gff_liftover)
subprocess.check_output(remap_command, shell=True)

gff_liftover = file(gff_liftover, 'r')

# Replace original coordinates
if vcf == True:
    gff = file(gff_temp, 'r')
    proc = Popen("bcftools view %s" % sys.argv[1], stdout=PIPE, stdin=PIPE, shell=True)
    for line in proc.stdout:
        line = line.replace("\n", "")
        if line.startswith("#"):
            pipe_out(line)
        else:
            # Add checks
            l = gff_liftover.readline().split("\t")
            pos_orig = l[9]
            pos_new = l[3]
            line = line.split("\t")
            if line[1] != pos_orig:
                raise Exception("Coordinates Off")
            else:
                line[1] = pos_new
                pipe_out('\t'.join(line))








