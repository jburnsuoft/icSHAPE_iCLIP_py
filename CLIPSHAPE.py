import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- Define helper functions ---
# This function reads an iCLIP .bed file and assigns proper column names.
def get_clip_coord(clipfile):
    clip_coord = pd.read_csv(clipfile, sep="\t", header=None)
    clip_coord.columns = ["chrom", "start", "end", "score1", "score2", "strand"]
    return clip_coord

# This function reads an icSHAPE coordinate .bed file and names columns.
def get_shape_coord(shape_file):
    shape_coord = pd.read_csv(shape_file, sep="\t", header=None)
    shape_coord.columns = ["chrom", "start", "end", "score"]
    return shape_coord

# This function checks each coordinate from iCLIP and finds overlaps with icSHAPE ranges.
def get_shape_range(zf_peaks, shape_coord, strand_sh, range_val):
    run_sites = pd.DataFrame()
    for x in range(len(zf_peaks)):
        has_run = True
        for i in range(-range_val, range_val + 1):
            # Calculate the site position offset by i.
            site = zf_peaks.iloc[x]["start"] + i
            # Check if the site exists with matching chromosome and strand.
            if not (site in shape_coord["start"].values) or \
               not (zf_peaks.iloc[x]["chrom"] in shape_coord[shape_coord["start"] == site]["chrom"].values) or \
               not (strand_sh == zf_peaks.iloc[x]["strand"]):
                has_run = False
                break
        if has_run:
            run_sites = run_sites.append(zf_peaks.iloc[x])
    return run_sites

# This function consolidates overlapping coordinates from both plus and minus strands.
def get_data_sites(clipfile, plus_shape, minus_shape, range_val):
    plus_sites = get_shape_range(clipfile, plus_shape, "+", range_val)
    minus_sites = get_shape_range(clipfile, minus_shape, "-", range_val)
    all_sites = pd.concat([plus_sites, minus_sites])
    return all_sites

# This function extracts and writes icSHAPE scores to output files based on provided coordinates.
def write_shape(pos_chrom, clus_pos, shape_data, init_range, run_size, seq_name, shape_folder_name):
    # Construct output file name
    seq_name = f"{shape_folder_name}{seq_name}.SHAPE"
    # Filter shape_data based on coordinate ranges
    shape_trimmed = shape_data[(shape_data["start"] >= (clus_pos - init_range)) & 
                               (shape_data["start"] <= (clus_pos + run_size)) & 
                               (shape_data["chrom"] == pos_chrom)]
    scores = shape_trimmed["score"]
    shape_out = pd.DataFrame({"index": range(1, len(scores) + 1), "scores": scores})
    # Write output using tab separator
    shape_out.to_csv(seq_name, sep="\t", index=False, header=False, quoting=False)

# --- Main function to run CLIPSHAPE ---
# This function combines iCLIP and icSHAPE data to identify binding sites, extract sequences,
# write individual .seq files, and finally generate an output fasta file.
def CLIPSHAPE(clipfile_name, shape_range, cluster_distance, plus_shape_name, minus_shape_name, out_fasta_name, seq_folder_name, shape_folder_name, ct_folder_name, db_folder_name):
    # Load reference genome in fasta format into a dictionary.
    chr_seq = SeqIO.to_dict(SeqIO.parse("GRCh37.p13.genome.fa", "fasta"))

    # Create output directories if they do not exist.
    os.makedirs(seq_folder_name[:-1], exist_ok=True)
    os.makedirs(shape_folder_name[:-1], exist_ok=True)
    os.makedirs(ct_folder_name[:-1], exist_ok=True)
    os.makedirs(db_folder_name[:-1], exist_ok=True)

    # Load iCLIP and icSHAPE coordinate files into DataFrame structures.
    full_clip = get_clip_coord(clipfile_name)
    shape_plus = get_shape_coord(plus_shape_name)
    shape_minus = get_shape_coord(minus_shape_name)

    # Identify overlapping coordinates with icSHAPE scores.
    clip = get_data_sites(full_clip, shape_plus, shape_minus, shape_range)
    if clip.empty:
        return "No binding sites with shape data found"

    # Sort and calculate differences between subsequent coordinates.
    clip = clip.sort_values(by=["chrom", "strand", "start"])
    clip["diff"] = clip.groupby(["chrom", "strand"])["start"].diff().shift(-1).fillna(0)

    # Prepare lists for sequences and their names.
    seq_list = []
    name_vec = []
    n = 0  # Used for skipping indices based on clustering.
    iter = 0

    # Loop over coordinates in the consolidated clip DataFrame.
    for i in range(len(clip)):
        if n > 0:
            n -= 1
            continue
        chrom = clip.iloc[i]["chrom"]
        # Determine the reference key for chromosome lookup.
        if chrom != "chrM":
            chrom_fasta = f"{chrom} {chrom[3:5]}"
        else:
            chrom_fasta = "chrM MT"
        strand = clip.iloc[i]["strand"]
        pos = clip.iloc[i]["start"]
        sequence = []  # Initialize sequence list for current binding site.
        name_vec.append(f"{chrom} {strand} {pos}")
        n = 0
        num = i
        count = shape_range
        # Group close coordinates based on cluster distance.
        while (num < len(clip)) and (not pd.isna(clip.iloc[num]["diff"])) and (clip.iloc[num]["diff"] < cluster_distance):
            count += shape_range - clip.iloc[num]["diff"]
            n += 1
            num = i + n
            name_vec[-1] = f"{name_vec[-1]} {clip.iloc[num]['start']}"
        iter += 1
        # Build the sequence across a sliding window defined by shape_range and count.
        x = -shape_range
        while x <= count:
            file_name_part = f"{chrom}.{pos}.{strand}"
            if strand == "+":
                sequence.append(str(chr_seq[chrom_fasta].seq[pos + x].complement()))
                if x == count:
                    # Write the icSHAPE score file when the window reaches its edge.
                    write_shape(chrom, pos, shape_plus, shape_range, count, file_name_part, shape_folder_name)
            else:
                sequence.append(str(chr_seq[chrom_fasta].seq[pos + x]))
                if x == count:
                    write_shape(chrom, pos, shape_minus, shape_range, count, file_name_part, shape_folder_name)
            # Convert any DNA to RNA nucleotide convention.
            sequence = [s.replace("T", "U") for s in sequence]
            x += 1
        # Write the extracted sequence to a .seq file if within a valid coordinate range.
        if (pos - shape_range) > 1:
            if strand == "+":
                seq_list.append(sequence[::-1])
            else:
                seq_list.append(sequence)
            with open(f"{seq_folder_name}{file_name_part}.seq", "w") as f:
                f.write("".join(sequence))
    
    # Build a FASTA file containing all extracted sequences.
    records = [SeqRecord(Seq("".join(seq)), id=name, description="") for seq, name in zip(seq_list, name_vec)]
    SeqIO.write(records, out_fasta_name, "fasta")

# --- Main wrapper to allow execution from the command line ---
if __name__ == '__main__':
    import sys
    # Ensure all required parameters are provided.
    if len(sys.argv) != 11:
        sys.exit("Usage: python CLIPSHAPE <clipfile_name> <shape_range> <cluster_distance> <plus_shape_name> <minus_shape_name> <out_fasta_name> <seq_folder_name> <shape_folder_name> <ct_folder_name> <db_folder_name>")
    _, clipfile_name, shape_range, cluster_distance, plus_shape_name, minus_shape_name, out_fasta_name, seq_folder_name, shape_folder_name, ct_folder_name, db_folder_name = sys.argv
    # Convert numeric parameters to integers.
    shape_range = int(shape_range)
    cluster_distance = int(cluster_distance)
    # Run the main CLIPSHAPE function.
    CLIPSHAPE(clipfile_name, shape_range, cluster_distance, plus_shape_name, minus_shape_name, out_fasta_name, seq_folder_name, shape_folder_name, ct_folder_name, db_folder_name)