import subprocess
from create_interm_RS_files import *
from create_interm_ENS_files import *


jar_file = f'{curwd}/snpEff_v4_3t_core/snpEff/snpEff.jar'

# Function to perform snpEff database build that is necessary for VCF annotation
def run_snpEff_build(files, config_path, species):
    # Creating the directory that snpEff will look for the data files in
    data_folder = f'{curwd}/data/{species}'
    try:
        os.makedirs(data_folder)
        print(f"Directory '{data_folder}' created successfully.\n")
    except FileExistsError:
        print(f"Directory '{data_folder}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{data_folder}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

        # snpEff does best with absolute paths, so ensuring the symlinks point back to the absolute file path
    for file in files:
        if not os.path.isabs(file):
            files.remove(file)
            files.append(os.path.abspath(file))

    # A list of symbolic links will be used in the snpEff data directory to refer to the files generated in the
    # previous modules, instead of the actual files, for simplicity
    links = [f'{data_folder}/genes.gff',
             f'{data_folder}/protein.fa',
             f'{data_folder}/cds.fa',
             f'{data_folder}/sequences.fa']

    symlink_dict = dict(zip(files, links))

    for original, link in symlink_dict.items():
        file_shuffle = ['ln', '-s', original, link]
        subprocess.run(file_shuffle, capture_output=True, check=True)
        print(f'File linked successfully.')

        # List of strings that are used in Linux to run the snpEff build command, which will be concatenated in the
        # subprocess.run() command to execute the build

    snpeff_build_command = f'java -jar {jar_file} build -gff3 -v {species} -c {curwd}/{config_path} -d'
    subprocess.run(snpeff_build_command, capture_output=True, shell = True, check=True)

    return print(f'\nsnpEff build command executed successfully.')


def run_snpEff_ann(VCF, species, config_path):
    vcf_path = os.path.abspath(VCF)
    annotated_file_name = f'{curwd}/ANN_{os.path.basename(VCF)}'

    snpeff_ann_command = f'java -Xmx8g -jar {jar_file} -c {curwd}/{config_path} -v {species} {vcf_path} > {annotated_file_name}'
    subprocess.run(snpeff_ann_command, capture_output=True, check=True, shell=True)

    print(f'\nAnnotation executed successfully. Your annotated VCF can be found at {annotated_file_name}.')
    return None