import argparse
import process_gff as prs
from process_gff import *
import create_interm_RS_files as RS
import create_interm_ENS_files as EN
import create_interm_JGI_files as JGI
import snpEff_exec


cl_parser = argparse.ArgumentParser(prog = 'Broom',
                                    description = 'A genomic data file preprocessing tool')
cl_parser.add_argument("--gff", required = True, help = "Path to the GFF file")
cl_parser.add_argument("--genome_fasta", required = True, help = "Path to the genome FASTA file")
cl_parser.add_argument("--repository", required = True, choices = ["Ensembl", "RefSeq","JGI-MP","JGI-Phy"],
                       help = "Data source. "
                              "JGI has two sub-repositories: Phytozome and Mycocosm/Phycocosm. This needs to"
                              "be specified as JGI-Phy or JGI-MP")
cl_parser.add_argument("--genome_name", required = True,
                       help = "Name of species; should be the same as the one listed in the config file")
cl_parser.add_argument("--vcf", required = False, help = "VCF you'd like to annotate with snpEff")
cl_parser.add_argument("-a", required = False, action='store_true',
                       help = "All-features flag, when turned on will result in all features being in reconstituted GFF,"
                              " regardless of if snpEff or any other tools will analyze them or not.")
cl_parser.add_argument("--excl", required = False, nargs='*', default = [],
                        help = "Input is the sequence/chromosome ID of any chromosomes you'd like to exclude from"
                               "analysis, i.e - mitochondria or chloroplast.")
cl_parser.add_argument("--config", required = False, help = "Path to the config file containing the "
                                                           "information necessary per the snpEff documentation. The genome "
                                                           "name you use in this file must be the same one passed for the"
                                                           "command line argument. Ensure any other edits particular to"
                                                           " your species are completed. Find relevant information at: "
                                                           "https://pcingola.github.io/SnpEff/snpeff/build_db/#step-1-configure-a-new-genome")


user_inputs = cl_parser.parse_args()

gff = user_inputs.gff
fasta_file = user_inputs.genome_fasta
repo = user_inputs.repository
species = user_inputs.genome_name
vcf = user_inputs.vcf
unwanted_ids = user_inputs.excl
config_path = user_inputs.config
all_features = user_inputs.a

if user_inputs.vcf and user_inputs.config is None:
    cl_parser.error(f"snpEff annotation requires a config file. Find instructions for how to format the config file"
                    f" at https://pcingola.github.io/SnpEff/snpeff/build_db/#step-1-configure-a-new-genome")

if user_inputs.config and user_inputs.vcf is None:
    cl_parser.error("snpEff annotation requires a VCF. Please provide one, along with a config file.")

if repo == "RefSeq":

    fileset = []

    process_gff(gff, species=species, repo=repo, unwanted_ids=unwanted_ids)

    RS.featureQC(accessors, repo=repo, species=species)

    unsorted_gff = RS.make_gff(accessors, species=species, repo=repo, gff=gff, all_features=all_features)

    fileset.append(RS.sort_gff(unsorted_gff, gff))

    fileset.extend(RS.make_fastas(fasta_file, accessors, repo=repo,species=species))

    fileset.append(fasta_file)

    if vcf:
        snpEff_exec.run_snpEff_ann(vcf, species=species, config_path=config_path)

if repo == "Ensembl":

    fileset = []

    process_gff(gff, species=species, repo=repo,unwanted_ids=unwanted_ids)

    EN.featureQC(accessors, repo=repo, species=species)

    unsorted_gff = EN.make_gff(accessors, species=species, repo=repo, gff=gff, all_features=all_features)

    fileset.append(EN.sort_gff(unsorted_gff, gff))

    fileset.extend(EN.make_fastas(fasta_file, accessors, repo=repo, species=species))

    fileset.append(fasta_file)

    if vcf:
        snpEff_exec.run_snpEff_ann(vcf, species=species, config_path=config_path)

if repo == "JGI-Phy":

    fileset = []

    process_gff(gff, species=species, repo=repo,unwanted_ids=unwanted_ids)

    JGI.featureQC_phy(accessors, repo=repo, species=species)

    unsorted_gff = JGI.make_gff_phy(accessors, species=species, repo=repo, gff=gff)

    fileset.append(JGI.sort_gff(unsorted_gff, gff))

    fileset.extend(JGI.make_fastas(fasta_file, accessors, repo=repo, species=species))

    fileset.append(fasta_file)

    if vcf:
        snpEff_exec.run_snpEff_ann(vcf, species=species, config_path=config_path)


if repo == "JGI-MP":

    fileset = []

    process_gff(gff, species=species, repo=repo,unwanted_ids=unwanted_ids)

    JGI.featureQC_mp(accessors, repo=repo, species=species)

    unsorted_gff = JGI.make_gff_mp(accessors, species=species, repo=repo, gff=gff)

    fileset.append(JGI.sort_gff(unsorted_gff, gff))

    fileset.extend(JGI.make_fastas(fasta_file, accessors, repo=repo, species=species))

    fileset.append(fasta_file)

    if vcf:
        snpEff_exec.run_snpEff_ann(vcf, species=species, config_path=config_path)