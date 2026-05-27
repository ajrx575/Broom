import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from process_gff import *
import warnings
from Bio import BiopythonWarning
import check_compression as cc


genetypes = ['gene', 'transposable_element_gene']
tscript_prot_types = ['CDS', 'exon', 'mRNA']
utrs = ['three_prime_UTR', 'five_prime_UTR']

curwd = os.getcwd()
deleted = set()

def featureQC_mp(db_dict, species, repo):
    with (open(f'{repo}_{species}_misformatted_gff_features.txt', 'a+') as gff_errors):
        for source, gu_db in db_dict.items():

            g_obj = list(gu_db.all_features())

            for feat in g_obj:

                if feat.id in deleted:
                    continue

                if feat.start > feat.end:
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the start position is greater")
                    gff_errors.write(
                        f" than the end position. Deleting feature from database, and excluding feature & its children")

                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)
                    continue

                if feat.strand not in ('+', '-', '.'):
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the strand is not +, -, or period.")
                    gff_errors.write(
                        f" Deleting feature from database and relatives from database, and excluding from minimized file."
                        f"\n{feat}\n")
                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)

                    continue

                if feat.frame not in ('0', '1', '2', '.'):
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the frame is not 0, 1, or 2.")
                    gff_errors.write(
                        f" Deleting feature and relatives from database, and from minimized file."
                        f"\n{feat}\n")

                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)

                    continue

                try:
                    if feat.featuretype in genetypes:
                        if feat.attributes['ID'][0]:
                            continue

                    elif feat.featuretype in tscript_prot_types:
                        if feat.featuretype != 'mRNA':
                            if feat.attributes['ID'][0] and feat.attributes['Parent'][0]:
                                continue

                        elif feat.attributes['ID'][0] and feat.attributes['Parent'][0] and feat.attributes['proteinId'][0]:
                            continue

                    elif feat.featuretype in utrs:
                        if feat.attributes['Parent'][0]:
                            continue

                except KeyError:
                    missing_keys = []

                    if 'ID' not in feat.attributes:
                        missing_keys.append('ID')

                        parentID = feat.attributes['Parent'][0]
                        all_parents = []
                        for f in g_obj:
                            if f.attributes['ID'][0] == parentID:
                                all_parents.append(f)
                        for item in all_parents:
                            deleted.add(item.id)
                            for child in gu_db.children(item.id):
                                deleted.add(child.id)

                            gu_db.delete(gu_db.children(item.id))
                            gu_db.delete(item)

                    if feat.featuretype in tscript_prot_types and 'Parent' not in feat.attributes:
                        missing_keys.append('Parent')

                        deleted.add(feat.id)
                        for child in gu_db.children(feat.id):
                            deleted.add(child.id)

                    if feat.featuretype == 'mRNA' and 'proteinId' not in feat.attributes:
                        missing_keys.append('ProteinID')

                        parentID = feat.attributes['Parent'][0]
                        all_parents = []
                        for f in g_obj:
                            if f.attributes['ID'][0] == parentID:
                                all_parents.append(f)
                        for item in all_parents:
                            deleted.add(item.id)
                            for child in gu_db.children(item.id):
                                deleted.add(child.id)

                            gu_db.delete(gu_db.children(item.id))
                            gu_db.delete(item)

                    gff_errors.write(f"#The following feature and its relatives will be excluded from the minimized ")
                    gff_errors.write(f"file for missing listed attribute(s):{', '.join(missing_keys)}:\n{feat}\n")

    print(f'\nFind erroneous features at {os.path.abspath(gff_errors.name)}')
    return deleted



def make_gff_mp(db_dict, species, repo, gff):

    if cc.comp_check(gff):
        handler = gzip.open(gff, 'rt', encoding='utf-8')
    else:
        handler = open(gff, 'r', encoding='utf-8')

    with (open(f'{repo}_{species}_interm.gff', 'a+') as file):

        for line in handler:
            if line.startswith('##gff'):
                file.write(line)
            if line.startswith('#NCBI'):
                file.write(line)
            if line.startswith('#created'):
                file.write(line)

        for source, gu_db in db_dict.items():
            g_obj = list(gu_db.all_features())

            for feat in g_obj:

                if feat.id in deleted:
                    continue

                if feat.featuretype in genetypes:
                    file.write(
                        f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                        f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]}\n")

                elif feat.featuretype in tscript_prot_types:
                    if feat.featuretype != 'mRNA':
                        file.write(
                            f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                            f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]};"
                            f"Parent={feat.attributes['Parent'][0]}\n")

                    else:
                        file.write(
                            f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                            f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]};"
                            f"Parent={feat.attributes['Parent'][0]};ProteinID={feat.attributes['proteinId'][0]}\n")

                elif feat.featuretype in utrs:
                    file.write(
                        f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                        f"{feat.strand}\t{feat.frame}\tParent={feat.attributes['Parent'][0]}\n")

    gff_path = os.path.abspath(file.name)

    print(f'Intermediate file creation completed for {species}.\n\n'
        f'Your GFF file will now be sorted.\n\n')
    return gff_path


def featureQC_phy(db_dict, species, repo,):
    with(open(f'{repo}_{species}_misformatted_gff_features.txt', 'a+') as gff_errors):
        for source, gu_db in db_dict.items():

            g_obj = list(gu_db.all_features())

            for feat in g_obj:

                if feat.id in deleted:
                    continue

                if feat.start > feat.end:
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the start position is greater")
                    gff_errors.write(
                        f" than the end position. Deleting feature from database, and excluding feature & its children")

                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)
                    continue

                if feat.strand not in ('+', '-', '.'):
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the strand is not +, -, or period.")
                    gff_errors.write(
                        f" Deleting feature from database and relatives from database, and excluding from minimized file."
                        f"\n{feat}\n")
                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)

                    continue

                if feat.frame not in ('0', '1', '2', '.'):
                    gff_errors.write(
                        f"#The following feature is not GFF format compliant because the frame is not 0, 1, or 2.")
                    gff_errors.write(
                        f" Deleting feature and relatives from database, and from minimized file."
                        f"\n{feat}\n")

                    parentID = feat.attributes['Parent'][0]
                    all_parents = []
                    for f in g_obj:
                        if f.attributes['ID'][0] == parentID:
                            all_parents.append(f)

                    for item in all_parents:
                        deleted.add(item.id)
                        for child in gu_db.children(item.id):
                            deleted.add(child.id)

                        gu_db.delete(gu_db.children(item.id))
                        gu_db.delete(item)

                    continue

                try:
                    if feat.featuretype in genetypes:
                        if feat.attributes['ID'][0]:
                            continue

                    elif feat.featuretype in tscript_prot_types:
                        if feat.featuretype != 'CDS':
                            if feat.attributes['ID'][0] and feat.attributes['Parent'][0]:
                                continue

                        elif feat.attributes['ID'][0] and feat.attributes['Parent'][0] and feat.attributes['pacid'][0]:
                            continue

                    elif feat.featuretype in utrs:
                        if feat.attributes['Parent'][0]:
                            continue

                except KeyError:
                    missing_keys = []

                    if 'ID' not in feat.attributes:
                        missing_keys.append('ID')

                        parentID = feat.attributes['Parent'][0]
                        all_parents = []
                        for f in g_obj:
                            if f.attributes['ID'][0] == parentID:
                                all_parents.append(f)
                        for item in all_parents:
                            deleted.add(item.id)
                            for child in gu_db.children(item.id):
                                deleted.add(child.id)

                            gu_db.delete(gu_db.children(item.id))
                            gu_db.delete(item)

                    if feat.featuretype in tscript_prot_types and 'Parent' not in feat.attributes:
                        missing_keys.append('Parent')

                        deleted.add(feat.id)
                        for child in gu_db.children(feat.id):
                            deleted.add(child.id)

                    if feat.featuretype == 'CDS' and 'pacid' not in feat.attributes:
                        missing_keys.append('ProteinID')

                        parentID = feat.attributes['Parent'][0]
                        all_parents = []
                        for f in g_obj:
                            if f.attributes['ID'][0] == parentID:
                                all_parents.append(f)
                        for item in all_parents:
                            deleted.add(item.id)
                            for child in gu_db.children(item.id):
                                deleted.add(child.id)

                            gu_db.delete(gu_db.children(item.id))
                            gu_db.delete(item)

                    gff_errors.write(f"#The following feature and its relatives will be excluded from the minimized ")
                    gff_errors.write(f"file for missing listed attribute(s):{', '.join(missing_keys)}:\n{feat}\n")

    print(f'\nFind erroneous features at {os.path.abspath(gff_errors.name)}')
    return deleted



def make_gff_phy(db_dict, species, repo, gff):

    if cc.comp_check(gff):
        handler = gzip.open(gff, 'rt', encoding='utf-8')
    else:
        handler = open(gff, 'r', encoding='utf-8')

    with (open(f'{repo}_{species}_interm.gff', 'a+') as file):

        for line in handler:
            if line.startswith('##'):
                file.write(line)

        for source, gu_db in db_dict.items():
            g_obj = list(gu_db.all_features())

            for feat in g_obj:

                if feat.id in deleted:
                    continue

                if feat.featuretype in genetypes:
                    file.write(
                        f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                        f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]}\n")

                elif feat.featuretype in tscript_prot_types:
                    if feat.featuretype != 'CDS':
                        file.write(
                            f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                            f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]};"
                            f"Parent={feat.attributes['Parent'][0]}\n")

                    else:
                        file.write(
                            f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                            f"{feat.strand}\t{feat.frame}\tID={feat.attributes['ID'][0]};"
                            f"Parent={feat.attributes['Parent'][0]};ProteinID={feat.attributes['pacid'][0]}\n")

                elif feat.featuretype in utrs:
                    file.write(
                        f"{feat.seqid}\t.\t{feat.featuretype}\t{feat.start}\t{feat.end}\t.\t"
                        f"{feat.strand}\t{feat.frame}\tParent={feat.attributes['Parent'][0]}\n")

    gff_path = os.path.abspath(file.name)

    print(f'\nIntermediate file creation completed for {species}.\n\n'
          f'Your GFF file will now be sorted.\n\n')
    return gff_path


def sort_gff(gff_path, original_gff):
    header_lines = []
    col_names = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    matches = ['start', 'end', 'strand']

    with open(gff_path,'r') as file:
        for line in file:
            if line.startswith('#'):
                header_lines.append(line)

    ori_df = pd.read_csv(original_gff, sep='\t', names=col_names, comment='#', low_memory=False)
    int_df = pd.read_csv(gff_path, sep='\t', names=col_names, comment='#', low_memory=False)

    merged = pd.merge(ori_df, int_df, on=matches, how='left')
    merged = merged.dropna()
    merged = merged.iloc[:, [0, 1, 2, 5, 7, 8, 9, 10, 11, 3, 4, 12, 6, 13, 14]]
    reduced_data = merged.iloc[:, 6:]

    x = os.path.basename(gff_path).replace('_interm', '')
    sorted_gff_path = f'{curwd}/{x}'

    with open(sorted_gff_path, 'w') as final:
        final.writelines(header_lines)
        reduced_data.to_csv(final, sep='\t', index=False, header=False)
        final.write(f'###')

    print(f'Find your sorted GFF at {os.path.abspath(sorted_gff_path)}\n\n')

    return sorted_gff_path


def make_fastas(genome_fasta, db_dict, species, repo):

    if cc.comp_check(genome_fasta):
        handler = gzip.open(genome_fasta, 'rt', encoding='utf-8')
    else:
        handler = open(genome_fasta, 'r', encoding='utf-8')

    if repo == 'JGI-MP':
        genome_sequences = SeqIO.to_dict(SeqIO.parse(handler, "fasta-pearson"))
    else:
        genome_sequences = SeqIO.to_dict(SeqIO.parse(handler, "fasta"))

    protein_records = []
    cds_records = []

    for source, gu_db in db_dict.items():

        desiredfeats = list(gu_db.features_of_type("mRNA"))

        for tscript in desiredfeats:
            identifier = tscript.id
            cds_features = list(gu_db.children(identifier, featuretype="CDS", order_by='start'))

            if not cds_features:
                continue

            seq_parts = [genome_sequences[cds.seqid].seq[cds.start - 1:cds.end] for cds in cds_features]

            dna_seq = Seq('').join(seq_parts)

            if cds_features[0].strand == "-":
                dna_seq = dna_seq.reverse_complement()

            cds_record = SeqRecord(
                dna_seq,
                id=cds_features[0].attributes['Parent'][0],
                description=f"")

            cds_records.append(cds_record)

            with (warnings.catch_warnings(record=True) as w):
                warnings.simplefilter("always", BiopythonWarning)
                protein_seq = dna_seq.translate(to_stop=True)

            with open(f'{repo}_{species}_fasta_warnings.txt', 'a') as fasta_errors:
                for warning in w:
                    if "Partial codon" in str(warning.message):
                        fasta_errors.write(f" Warning: Partial codon in feature:\n"
                                           f"{cds_features[0]}\n")
                fasta_errors_path = os.path.abspath(fasta_errors.name)

            protein_record = SeqRecord(
                protein_seq,
                id=cds_features[0].attributes['Parent'][0],
                description=f"")
            protein_records.append(protein_record)

        # Create output files
    pfa_name = f"{repo}_{species}_protein.fa"
    cfa_name = f"{repo}_{species}_cds.fa"
    SeqIO.write(protein_records, pfa_name, "fasta")
    SeqIO.write(cds_records, cfa_name, "fasta")

    protein_path = os.path.abspath(pfa_name)
    cds_path = os.path.abspath(cfa_name)


    print(f'\nIntermediate FASTA files created.\n\n'
          f'Find erroneous features at {fasta_errors_path}\n\n'
          f'Protein FASTA can be found at {protein_path}\n\n'
          f'CDS FASTA can be found at {cds_path}\n\n')
    return protein_path, cds_path
