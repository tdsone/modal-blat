import modal

app = modal.App("razers")

image = (
    modal.Image.debian_slim()
    .apt_install("wget", "libcurl4")
    .run_commands(
        "wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat",
        "chmod +x blat",
    )
    .pip_install("pandas", "biopython")
)

blat_result_vol = modal.Volume.from_name("blat-results", create_if_missing=True)


def get_sequences_df(fasta_file: str, gff_file: str):
    import os

    print(os.listdir("/root"))

    import pandas as pd
    from Bio import SeqIO

    gff_data = pd.read_csv(
        gff_file,
        sep="\t",
        header=None,
        names=[
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )

    gff_data["id"] = gff_data.apply(
        lambda row: row["attributes"].split(";")[0].split("_")[0], axis=1
    )

    loss_regions = gff_data[gff_data["feature"] == "loss_region"]
    sample_regions = gff_data[gff_data["feature"] == "sample_region"]

    # Merge loss_region and sample_region on id
    merged_regions = pd.merge(
        loss_regions, sample_regions, on="id", suffixes=("_loss", "_sample")
    )

    samples = merged_regions[
        [
            "id",
            "chr_loss",
            "strand_loss",
            "start_loss",
            "end_loss",
            "start_sample",
            "end_sample",
        ]
    ]

    # Open fasta file and add the corresponding sequences to the dataframe
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    def extract_sequence(row, seq_dict, start_col, end_col, chr_col):
        seq_id = row[chr_col]
        start = row[start_col]
        end = row[end_col]
        return str(seq_dict[seq_id].seq[start - 1 : end])

    samples["loss_sequence"] = samples.apply(
        extract_sequence, args=(sequences, "start_loss", "end_loss", "chr_loss"), axis=1
    )

    samples["sample_sequence"] = samples.apply(
        extract_sequence,
        args=(sequences, "start_sample", "end_sample", "chr_loss"),
        axis=1,
    )

    return samples


@app.function(image=image, volumes={"/data": blat_result_vol}, timeout=600, cpu=16)
def run_razers():
    import os

    print("Fetching sequences...")

    seq_df = get_sequences_df(
        gff_file="/root/src/output_5000_2e-10.gff", fasta_file="/root/src/R64-1-1.fa"
    )

    print("Done fetching sequences.")
    print("Finding homologies...")

    os.makedirs("/root/src/fasta")
    os.makedirs("/data/blat", exist_ok=True)

    # Compare each sequence against the reference genome using razers
    for _, sample in seq_df.iterrows():
        # Write sample full_region sequence to fasta
        id = sample["id"].split("=")[1]
        sample_fasta_filename = f"/root/src/fasta/sample_{id}.fasta"
        with open(sample_fasta_filename, "w") as fasta_file:
            fasta_file.write(f">{id}\n{sample['sample_sequence']}\n")

        command = f"/blat /root/src/R64-1-1.fa {sample_fasta_filename} /data/blat/out_{id}.blast8 -out=blast8"

        # Execute command and print all logs and errors to stdout
        import subprocess

        process = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        for line in process.stdout:
            print(line.decode(), end="")
        process.wait()

    print("Done finding homologies")

    pass


@app.local_entrypoint()
def local():
    run_razers.remote()
