import gffutils
import check_compression as cc
import gzip
import os


accessors = {}

def process_gff(gff, species, repo, unwanted_ids):

    allfiles = set()

    if cc.comp_check(gff):
        handler = gzip.open(gff, 'rt', encoding='utf-8')
    else:
        handler = open(gff, 'r', encoding='utf-8')

    parts = []
    for line in handler:
        if not line.startswith('#'):
            parts.append(line.strip().split('\t')[1])

    handler.seek(0)
    unique_annot_pipelines = set(parts)

    for line in handler:
        if not line.startswith('#'):
            parts = line.strip().split('\t')
            if parts[0] in unwanted_ids:
                print(f'Skipping unwanted chromosome: {parts[0]}\n')
                continue
            for item in unique_annot_pipelines:
                if parts[0] not in unwanted_ids and parts[1] == item:
                    desired_annot_file = f'{repo}_{species}_{item}.gff'
                    allfiles.add(desired_annot_file)
                    with open(desired_annot_file, 'a') as i:
                        i.write(line)

    handler.close()

    allfiles = sorted(allfiles)
    print(f'\n\nSource-specific GFFs created.\n\n')

    db_names = []
    for newfile in allfiles:
        dbfn = f'{os.path.splitext(os.path.basename(newfile))[0]}.sqlite'
        db_names.append(dbfn)
        gffutils.create_db(newfile, dbfn=dbfn,
                           force = True,
                           keep_order=True,
                           merge_strategy='create_unique',
                           force_gff = True )
        print('gffutils DB created.')

    db_dict = dict(zip(sorted(db_names), sorted(unique_annot_pipelines)))
    for key, val in db_dict.items():
        accessors[val] = gffutils.FeatureDB(key, keep_order=True)

    return accessors