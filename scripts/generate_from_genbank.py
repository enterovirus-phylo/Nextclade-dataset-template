#!/bin/env python3
import os
import sys
from Bio import SeqIO
from collections import defaultdict
import re

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--reference', required=True, type=str, help='Genbank accession of reference sequence (will be fetched from genbank)')
    parser.add_argument('--output-dir', required=True, type=str, help='Output directory')
    return parser.parse_args()

def get_entrez_email():
    """
    NCBI requires a contact address with every Entrez request and throttles or
    blocks anonymous traffic. Read it from the EMAIL environment variable, or
    from a .env file in the project root (which .gitignore keeps out of git).
    """
    try:
        from dotenv import load_dotenv, find_dotenv
        load_dotenv(find_dotenv())
    except ImportError:
        # python-dotenv is optional; EMAIL may equally be exported by hand.
        pass

    email = os.environ.get("EMAIL", "").strip()
    return email or "you@example.org"


def get_reference_sequence(accession):
    from Bio import Entrez

    Entrez.email = get_entrez_email()
    print(f"Fetching reference sequence {accession} from GenBank")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    try:
        return SeqIO.read(handle, "genbank")
    except ValueError:
        sys.exit(
            f"\nCould not read a GenBank record for '{accession}'.\n"
            "Check that the accession exists and includes a version suffix, e.g. AY426531.1\n"
        )

def get_gff(accession):
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={accession}"
    import urllib
    return [x.decode() for x in urllib.request.urlopen(url).readlines()]

def annotate_sequence(seq, cds):
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation, SimpleLocation, CompoundLocation
    strands = {'+':1, '-':-1}
    new_features = [feat for feat in seq.features if feat.type=='source']
    for cdsid, segments in cds.items():
        if len(segments)==0:
            continue
        elif len(segments)==1:
            loc = FeatureLocation(int(segments[0][0][3])-1, int(segments[0][0][4]),
                                  strand = strands.get(segments[0][0][6], 0))
        else:
            loc = CompoundLocation([SimpleLocation(int(seg[0][3])-1, int(seg[0][4]),
                                  strand = strands.get(seg[0][6], 0)) for seg in segments])

        feat = SeqFeature(location = loc, type=segments[0][0][2],
                          qualifiers=segments[0][1], id=cdsid)
        # to appease augur, overwrite 'locus_tag' with the 'Name' attribute used in the GFF file to identify the CDS
        feat.qualifiers['locus_tag'] = feat.qualifiers['Name']
        new_features.append(feat)
    seq.features = new_features
    return seq

FORBIDDEN_NAME_CHARS = ' ,.;:()[]'


def clean_qualifier(value):
    """
    Strip the boilerplate that GenBank product names tend to carry, so that
    "VP1 protein" and "3C protease (3C)" come out as "VP1" and "3C".
    Shorter names make for readable amino acid mutation labels.
    """
    # Require a separator around the boilerplate word, so that "VP1 protein"
    # loses its suffix but "polyprotein" survives intact.
    value = re.sub(r'^(protein|protease|protien)[\s_|]+', '', value, flags=re.IGNORECASE)
    value = re.sub(r'[\s_|]+(protein|protease|protien)$', '', value, flags=re.IGNORECASE)
    value = re.sub(r'\s*\([\w\d]+\)\s*$', '', value)
    value = re.sub(r'__\d[A-D]_.*', '', value)
    return value.strip()


def check_name(new_name):
    """
    Normalise a CDS name: drop the boilerplate, then replace any character
    Nextclade and augur cannot handle in a feature name.
    """
    new_name = clean_qualifier(new_name)

    if any(c in new_name for c in FORBIDDEN_NAME_CHARS):
        print(f"WARNING: the CDS name '{new_name}' contains characters that are not allowed "
              f"in a feature name ('{FORBIDDEN_NAME_CHARS}').")
        for c in FORBIDDEN_NAME_CHARS:
            new_name = new_name.replace(c, '_')
        # Collapse the runs of underscores that e.g. ", " leaves behind.
        new_name = re.sub(r'_+', '_', new_name).strip('_')
        print(f"         It has been rewritten as '{new_name}'. Consider choosing a shorter name yourself.")

    if len(new_name) > 20:
        print(f"WARNING: the CDS name '{new_name}' is long; mutation labels will be unwieldy.")

    return new_name

if __name__=="__main__":
    args = parse_args()
    # Tolerate a trailing slash, so paths in the messages below read cleanly.
    args.output_dir = args.output_dir.rstrip("/") or "."
    reference = get_reference_sequence(args.reference)
    gff = get_gff(args.reference)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # write the reference sequence to the output directory
    SeqIO.write(reference, f"{args.output_dir}/reference.fasta", "fasta")

    # collect all cds and alternative annotations of protein sequences
    all_cds = defaultdict(lambda: defaultdict(list))
    for line in gff:
        if line[0]=='#':
            continue
        entries = line.strip().split('\t')
        if len(entries)<9:
            continue # invalid line

        # parse attributes and feature type and ID
        attributes = {x.split('=')[0]:x.split('=')[1] for x in entries[-1].split(';')}
        feature_type = entries[2]
        # IDs look like this: ID=id-NP_057850.1:133..363 where the part after to colon is the range in the translated sequence (exists in cases of mature protein annotations)
        feature_id = attributes['ID'].split(':')[0].split('-')[-1]
        if feature_type=='CDS':
            all_cds[feature_id]['CDS'].append([entries[:-1], attributes])
        elif feature_type in ['mature_protein_region_of_CDS', 'mat_peptide', 'mat_protein']:
            all_cds[feature_id]['mature_protein'].append([entries[:-1], attributes])


    annotation_choice = 0
    if any([len(feat)>1 for feat in all_cds.values()]):
        print("\nFor polypeptides annotated in multiple ways, you need to choose between different annotations (CDS or mature peptides/proteins).\n")
        while True:
            annotation_choice = input("For features with two annotations, use always CDS [1], always mature_protein [2], or pick on a case-by-case basis [0]: ")
            try:
                annotation_choice = int(annotation_choice)
                assert annotation_choice<3
                break
            except:
                print(f"You input '{annotation_choice}' is not valid. Please enter one of [0, 1, 2]")

    available_attributes_by_type = {'CDS':set(), 'mature_protein': set()}
    for annot, annot_set in available_attributes_by_type.items():
        for feat in all_cds.values():
            if annot in feat:
                for seg in feat[annot]:
                    if len(annot_set):
                        annot_set = annot_set.intersection(list(seg[1].keys()))
                    else:
                        annot_set = annot_set.union(list(seg[1].keys()))
        available_attributes_by_type[annot] = annot_set
        print(f"\nCommon attributes of all available {annot} annotations are:\n", annot_set, '\n')

    available_attributes = available_attributes_by_type['CDS'] if annotation_choice==1 else \
                           available_attributes_by_type['CDS'].intersection(available_attributes_by_type['mature_protein'])

    print("Note that CDS names should be short, unique, and recognizable as they are used to label amino acid mutations. Space, commas, colons, etc are not allowed.\n")
    while True:
        name_choice = input("Specify field to use as CDS name or leave empty for manual choice: ")
        if name_choice in available_attributes or name_choice=='':
            break
        else:
            print("invalid choice: ", name_choice)
            print("pick one of", available_attributes)

    streamlined_cds = {}
    names_by_id = {}
    name_count = defaultdict(int)
    # loop through all CDS and ask the user to pick one of the annotations
    # allow renaming of the CDS to user friendly names
    for cds_id, cds_sets in all_cds.items():
        cds_set = list(cds_sets.keys())[0]
        if annotation_choice==0:
            if len(cds_sets)>1:
                print(f"\nCDS {cds_id} is annotated in multiple ways:")
                for i, (cds_set, segments) in enumerate(cds_sets.items()):
                    print(f"\t{cds_set} with a total of {len(segments)} items [{i+1}]")
                choice = int(input("Please pick number in brackets to choose one (0 for omission): "))
                if choice:
                    cds_set = list(cds_sets.keys())[choice-1]
                else:
                    continue
        else:
            if len(cds_sets)>1:
                cds_set = list(cds_sets.keys())[annotation_choice-1]

        segments = cds_sets[cds_set]

        for segment in segments:
            segment_id = segment[1]["ID"]
            # only needed for one segment per CDS, skip if already in streamlined_cds
            if segment_id not in streamlined_cds:
                if name_choice:
                    new_name = segment[1][name_choice]
                else:
                    print(f"Attributes of the segment with ID='{segment_id}' are:")
                    for k,v in segment[1].items():
                        print(f'\t\t{k}:\t{v}')
                    new_name = input("Enter desired name (leave empty to drop): ")
                    if new_name=='': continue

                streamlined_cds[segment_id] = []
                new_name = check_name(new_name)
                name_count[new_name] += 1
                if name_count[new_name]>1:
                    print(f"\nWARNING: CDS name {new_name} is not unique, attaching index to make unique.\n")
                    new_name = f"{new_name}_{name_count[new_name]:03d}"

                names_by_id[segment_id] = new_name

            # if renamed and selected, add the segment to the list of segments
            if segment_id in names_by_id:
                new_entries, new_attributes = list(segment[0]), dict(segment[1])
                new_entries[2]='CDS'
                new_attributes['Name']=names_by_id[segment_id]
                if 'product' in new_attributes:
                    new_attributes['product'] = check_name(new_attributes['product'])
                if 'locus_tag' in new_attributes:
                    new_attributes['locus_tag'] = check_name(new_attributes['locus_tag'])
                if "Parent" in new_attributes: new_attributes.pop("Parent")
                streamlined_cds[segment_id].append([new_entries, new_attributes])

    gff_fname = f"{args.output_dir}/genome_annotation.gff3"
    print(f"\nWriting annotation to GFF file: {gff_fname}\n")
    # write the gff file as a simple text file line by line
    with open(gff_fname, "w") as f:
        # write the header and the region line
        for line in gff:
            entries = line.split('\t')
            if entries[0][0]=='#' or (len(entries)>8 and entries[2]=='region'):
                f.write(line)

        # write the CDS lines
        for cds in streamlined_cds:
            if len(streamlined_cds[cds])==0: continue
            print(f"Exporting CDS '{names_by_id[cds]}' with {len(streamlined_cds[cds])} segments.")
            for segment in streamlined_cds[cds]:
                attributes = ';'.join([f"{k}={v}" for k,v in sorted(segment[1].items(), key=lambda x:(x[0]!='Name', len(x[1])))])
                f.write('\t'.join(segment[0])+'\t' + attributes + '\n')

    gb_fname = f"{args.output_dir}/reference.gbk"
    print("\nWriting annotation to genbank file: ", gb_fname)
    reannotated_seq = annotate_sequence(reference, streamlined_cds)
    SeqIO.write(reannotated_seq, gb_fname, "genbank")

    print(
        "\nDone. Two files are ready for the Nextclade dataset:\n"
        f"  {args.output_dir}/reference.fasta\n"
        f"  {gff_fname}\n"
        "\nThe workflow also needs the annotated GenBank file in resources/, where\n"
        "augur ancestral looks for it:\n"
        f"  cp {gb_fname} resources/reference.gbk\n"
        "\nThen make sure the GENES list in the Snakefile matches the CDS names above."
    )