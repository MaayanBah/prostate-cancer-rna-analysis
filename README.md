# Prostate Cancer RNA Analysis

## Background

Prostate cancer arises in the prostate gland, a small walnut-shaped organ in males responsible for producing seminal
fluid vital for nourishing and transporting sperm. Though it's one of the most prevalent cancer types, many cases
progress slowly, remaining localized within the prostate without posing significant harm. However, certain aggressive
forms can rapidly metastasize.

Early detection is key, as prostate cancer confined to the gland offers the highest treatment success rates.
Our objective involves analyzing RNA-seq data from 14 primary prostate cancers and their corresponding normal
tissues from individuals of Chinese descent. Through this analysis, we aim to gain insights into the molecular
mechanisms underlying prostate cancer progression in this population.

## Step 1

Using fastq-dump, we downloaded 3,000,000 reads from each RNA-seq sample described in the table to the class
machine:

| Run       | SampleName | Condition |
|-----------|------------|-----------|
| ERR031027 | 1          | Normal    |
| ERR031028 | 1          | Tumor     |
| ERR031029 | 2          | Normal    |
| ERR031030 | 2          | Tumor     |
| ERR031031 | 3          | Normal    |
| ERR031032 | 3          | Tumor     |
| ERR031033 | 4          | Normal    |
| ERR299297 | 4          | Tumor     |
| ERR031035 | 5          | Normal    |
| ERR299298 | 5          | Tumor     |
| ERR299299 | 6          | Normal    |
| ERR031038 | 6          | Tumor     |
| ERR031039 | 7          | Normal    |
| ERR031040 | 7          | Tumor     |
| ERR031041 | 8          | Normal    |
| ERR031042 | 8          | Tumor     |
| ERR031043 | 9          | Normal    |
| ERR031044 | 9          | Tumor     |
| ERR031017 | 10         | Normal    |
| ERR031018 | 10         | Tumor     |
| ERR031019 | 11         | Normal    |
| ERR299295 | 11         | Tumor     |
| ERR299296 | 12         | Normal    |
| ERR031022 | 12         | Tumor     |
| ERR031023 | 13         | Normal    |
| ERR031024 | 13         | Tumor     |
| ERR031025 | 14         | Normal    |
| ERR031026 | 14         | Tumor     |

```bash
for file in ERR031027 ERR031028 ERR031029 ERR031030 ERR031031 ERR031032 ERR031033 ERR299297 ERR031035 ERR299298 \
ERR299299 ERR031038 ERR031039 ERR031040 ERR031041 ERR031042 ERR031043 ERR031044 ERR031017 ERR031018 ERR031019 \
ERR299295 ERR299296 ERR031022 ERR031023 ERR031024 ERR031025 ERR031026; do
    fastq-dump $file --split-files --gzip -X 3000000 &
done

# Wait for all background processes to finish
wait
```

I used the fastq-dump command on each run accession in the table above to download 3,000,000 reads of its data,
utilizing the -X parameter. Additionally, I employed --gzip to conserve disk space and --split-files for paired-end
data. The & option enabled the commands to run in parallel.

```bash
# Make a new directory for the files.
mkdir prostate_cancer_samples
mv ERR* prostate_cancer_samples
```

---

## Step 2

Cleaning the data is a must in this assignment. It's only natural – the data is over 10 years old and considered
ancient in computational biology time scales. Therefore, the question is not whether to clean it but how to do it.
We'll start by assessing the raw data quality using FastQC and MultiQC, followed by cleaning it with fastp.
Subsequently, we'll reevaluate the data's quality using FastQC and MultiQC to ensure it meets our standards.

```bash
# Quality check using fastqc 
for file in *.fastq.gz; do
    fastqc "$file" &
done
```

We generated a FastQC report for each file, and now we'll proceed to create a MultiQC report.

```bash
mkdir first_quality_reports
mv *.html /home/maayanbah/hw_5/first_quality_reports
mv *.zip /home/maayanbah/hw_5/first_quality_reports

# quality check using multiqc
multiqc /home/maayanbah/hw_5/first_quality_reports
```

I copied the files to my machine to view them. To do this, I used the following command:

```bash
exit
scp -r maayanbah@10.0.32.173:/home/maayanbah/hw_5/first_quality_reports first_quality_reports
```

After I copied the files, I used WinSCP to transfer them to my computer.
First, we'll review the fastqc reports. You can view screenshots from the report in the attached first_quality_reports
directory.

Basic Statistics: I checked that all the samples have fastq reports.
Sequence Quality Histograms: There are several things that I've noticed:

1. Two samples had low quality at the first 3-5 bases. Therefore, I decided to trim the first 5 bases.
2. ERR031025 had low quality in the middle of the file, since it appeared only in one file and in the middle of the
   sequence, and due to the large number of samples, I decided to remove this sample and its coupled sample. 
   Therefore, I removed ERR031025 and ERR031026.
3. At the last 24 bases, the average quality is lower than 25, I decided to cut these bases from all samples.

Per Base Sequence Content: we have a problematic result in the first 10 bp. it's not a good result but on its own,
it might not be enough for us to decide that we want to cut it. since the "Basic Statistics" results support the claim
that the first 5 bp are problematic, I've decided to trim it.
I also noticed that about 15 files had bad results throughout the file. I opened the fastqc reports and looked at the
"Per base sequence content" table since it's easier to look at. As you can see in the first_quality_reports directory,
we have more A and T than G and C. Although we expect to see about 25% of 
each nucleotide - sometimes the A/T content will be higher since G/C are rarer.
For this reason, I decided to keep these files and use them.

Per Sequence GC Content: The distribution of GC content is not very good in some of the files.

Per Base N Content: The ERR031025 file had poor results (it had a high N content). since I had already decided to remove
this sample, I ignored it. The other samples had excellent results.

Sequence Duplication Levels: The level of duplication is a bit high.

After reviewing the results, I decided to use Fastq in order to clean the samples.

```bash
# Function to run fastp command for each pair of files
run_fastp() {
    # This assigns the value of the first positional parameter ($1) passed to the function to the variable.
    local file="$1"
    # We remove the extension of the file.
    local base_name="${file%_1.fastq.gz}"
    
    # Running fastq
    fastp \
    --in1 "${base_name}_1.fastq.gz" \
    --in2 "${base_name}_2.fastq.gz" \
    --out1 "${base_name}_1.fastp.fastq.gz" \
    --out2 "${base_name}_2.fastp.fastq.gz" \
    --unpaired1 "${base_name}_1.unpaired.fastq" \
    --unpaired2 "${base_name}_2.unpaired.fastq" \
    --average_qual 30 \
    --trim_front1 5 \
    --trim_tail1 24 \
    --length_required 61
}

# Iterate over each pair of FASTQ files and run fastp in parallel
for file in *_1.fastq.gz; do
    run_fastp "$file" &
done

# Wait for all background processes to finish
wait

# Remove the unwanted samples: ERR031025 and ERR031026.
rm ERR031025_1.fastp.fastq.gz
rm ERR031025_2.fastp.fastq.gz
rm ERR031025_1.unpaired.fastq
rm ERR031025_2.unpaired.fastq
rm ERR031026_1.fastp.fastq.gz
rm ERR031026_2.fastp.fastq.gz
rm ERR031026_1.unpaired.fastq
rm ERR031026_2.unpaired.fastq
```

fastq clean options explained:

1. -in1: first input.
2. -in2: second input.
3. --out1: first output.
4. --out2: second output.
5. --unpaired1: an output file for unpaired reads (orphan reads).
6. --unpaired2: an output file for unpaired reads (for the second file).
7. --average_qual: minimum average base quality per read, I picked 30.
8. --trim_front1: number of bases to trim from the start of each read of the first input FastQ - I chose to cut 10 bp
   from the start in order to get rid of the 10 problematic nucleotides at the beginning (as seen in the
   "Per Base Sequence Content" table)
9. --trim_tail1: number of bases to trim from each end - I chose 24 since the quality decrease over time,
   so I trim the very end in order to get a higher quality.
10. --length_required: I chose to set it to be 56 in order to keep a uniform length distribution.


The files ending with _1.fastp.fastq.gz and _2.fastp.fastq.gz are the processed paired-end FASTQ files containing
the reads after trimming and quality filtering performed by fastp, I moved it to a new directory:
Files ending with _1.unpaired.fastq and _2.unpaired.fastq are the unpaired reads that were not merged into pairs
during the processing. This can happen if one of the reads in a pair did not meet the filtering criteria or if the
pair was unpaired due to adapter contamination or other issues.

```bash
mkdir clean_samples
mv *_1.fastp.fastq.gz /home/maayanbah/hw_5/clean_samples
mv *_2.fastp.fastq.gz /home/maayanbah/hw_5/clean_samples
mv *_1.unpaired.fastq /home/maayanbah/hw_5/clean_samples
mv *_2.unpaired.fastq /home/maayanbah/hw_5/clean_samples
```

We have 103 files as expected. now we will create fastqc reports again:

```bash
# Second quality check using fastqc 
for file in *.fastp.fastq.gz; do
    fastqc "$file" &
done
```

We generated a second FastQC report for each file, and now we'll proceed to create a MultiQC report.

```bash
mkdir second_quality_reports
mv *.html /home/maayanbah/hw_5/second_quality_reports
mv *.zip /home/maayanbah/hw_5/second_quality_reports

# quality check using multiqc
multiqc /home/maayanbah/hw_5/second_quality_reports
```

Once again I copied the files to my machine to view them using the following command:

```bash
exit
scp -r maayanbah@10.0.32.173:/home/maayanbah/hw_5/second_quality_reports second_quality_reports
```

After I copied the files, I used WinSCP to transfer them to my computer.
Now, we can review the fastqc reports. You can view the report in the second_quality_reports directory.


Basic Statistics: All files are present.

Sequence Quality Histograms: We got great results, the average quality is higher than 30.

Per Base Sequence Content: The results have improved, but they're not perfect. The beginning still doesn't show very
good results. However, if the other graphs are satisfactory, I'll overlook this issue. I chose not to remove the files
with high content of T and A because, as explained earlier, it's normal to observe higher A/T content in certain cases.

Per Sequence GC Content: The results are similar to before. It's not perfect, but good enough for the analysis.

Per Base N Content: excellent results.

Sequence Duplication Levels: The level of duplication is still moderately high.

These results are good enough for the analysis.

---

## Step 3

We'll use bwa mem to align the data to the human genome and save the alignment in BAM format, sorted, and indexed.
bwa mem is not the best tool for the job, but for the sake of simplicity, it will suffice.

Before aligning my reads I need to index the reference genome using bwa index:

```bash
# activate our conda environment
conda activate compgen

# download the human genome
mkdir -p ~/GenomicPractice/genomes/human/bwa_index
cd ~/GenomicPractice/genomes/human/bwa_index 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# prepare an index for the human genome
bwa index hg38.fa.gz
```

Now we'll align the reads, we'll limit the use of resources since it's a heavy command and we don't want the jobs
to be killed.

```bash
# Aligning the Reads
for file in *_1.fastp.fastq.gz; do    
    base_name=$(basename "$file" _1.fastp.fastq.gz)
    bwa mem /home/maayanbah/GenomicPractice/genomes/human/bwa_index/hg38.fa.gz \
    "${base_name}_1.fastp.fastq.gz" "${base_name}_2.fastp.fastq.gz" -t 8 > "${base_name}.sam"
    
done
```

We use samtools view to retain only properly-paired aligned reads (filtering for properly paired reads can improve the
quality of out results), and save the alignment in BAM format, sorted,
and indexed:

```bash
run_samtools() {
   local file="$1"
   # We remove the extension of the file.
   base_name=$(basename "$file" _1.fastp.fastq.gz)
   # Create a bam format out of the SAM file and also use the -F flag to keep only the reads that are properly paired
   # and the -b flag to get output in  BAM format
   samtools view -b -f 2 "${base_name}.sam" > "${base_name}.bam"
   # Create a sorted BAM file.
   samtools sort "${base_name}.bam" > "${base_name}.sorted.bam"
   # Index the sorted BAM file.
   samtools index "${base_name}.sorted.bam"
}

for file in *_1.fastp.fastq.gz; do
  run_samtools "$file" &
done
wait

# Concatenate all BAM files into a single file
(for file in *.bam; do
  samtools flagstat "$file"; 
done) > flagstats_reports.bam.flagstats
```

After this part I checked the flagstats report:
![image](https://github.com/MaayanBah/prostate-cancer-rna-analysis/assets/84293984/16d692c7-43cb-4483-9013-2190ff447f6b)
The flagstats report appears excellent: every read is properly paired, all reads have been successfully mapped,
and the majority are classified as primary reads, meaning they are of high quality. Additionally, there are no
duplicates.

## Step 4

In this step we'll use bcftools to call variants in each sample and then merge the individual VCF files by groups,
that is, one merged file for normal tissues and one for tumor ones. We'll use only high-quality variants (e.g., with
good quality scores and sufficient coverage).

```bash
run_bcftools() {
   local file="$1"
   base_name=$(basename "${file}" _1.fastp.fastq.gz)
   # Call variants
   # The -OU flag specifies the output format as uncompressed BCF
   # The -f flag specifies the reference genome sequence file 
   # The -mv flag: -m switch tells the program to use the default calling method, the -v option asks to output only variant sites 
   # The -Ob flag: the -O option selects the output format. In this example we chosen binary compressed BCF.
   # The -o specifies the output file name.
   bcftools mpileup -Ou -f /home/maayanbah/GenomicPractice/genomes/human/bwa_index/hg38.fa \
   "${base_name}".sorted.bam | bcftools call -mv -Ob -o "${base_name}".bcf
   bcftools sort -o "${base_name}"_sorted.bcf "${base_name}".bcf
   bcftools index "${base_name}"_sorted.bcf
    # We want to choose reads with high quality:
    # DP: typically indicates the number of reads aligned to that position. higher depth values suggest more confidence
    # in the accuracy of the variant call.
    # QUAL: The QUAL field represents the probability that the variant is real.
   bcftools view -o "${base_name}"_s_filtered.vcf -i 'DP >= 10 && QUAL >= 30' "${base_name}"_sorted.bcf
   bgzip "${base_name}"_s_filtered.vcf
   bcftools index -f "${base_name}"_s_filtered.vcf.gz
}

for file in *_1.fastp.fastq.gz; do
  run_bcftools "${file}" &
done
wait

# Now we'll merge the files to two groups - tumor samples and normal samples
tumor_command="bcftools merge -o tumor_samples.vcf"
normal_command="bcftools merge -o normal_samples.vcf"

tumor_files=(
    "ERR031028"
    "ERR031030"
    "ERR031032"
    "ERR299297"
    "ERR299298"
    "ERR031038"
    "ERR031040"
    "ERR031042"
    "ERR031044"
    "ERR031018"
    "ERR299295"
    "ERR031022"
    "ERR031024"
)
normal_files=(
   "ERR031027"
   "ERR031029"
   "ERR031031"
   "ERR031033"
   "ERR031035"
   "ERR299299"
   "ERR031039"
   "ERR031041"
   "ERR031043"
   "ERR031017"
   "ERR031019"
   "ERR299296"
   "ERR031023"
)

# Add each tumor file to the command
for sample_name in "${tumor_files[@]}"; do
  tumor_command+=" ${sample_name}_s_filtered.vcf.gz"
done

# Add each sample file to the command
for sample_name in "${normal_files[@]}"; do
  normal_command+=" ${sample_name}_s_filtered.vcf.gz"
done

# Run the command
eval "${tumor_command}"
eval "${normal_command}"
```

**stopped here**

## Step 5

Now, we would like to filter and retain only somatic and novel SNP variants. This means we are interested in SNPs that
are present in tumor tissues but not in normal tissues (somatic variants) and SNPs that have not been previously
identified (novel variants).

We make the following assumptions:

1. Consider a position X at chromosome Y in any merged VCF file. If any alternative alleles in this position are SNPs,
   we consider the position a tissue variant.
2. We only compare SNPs by their coordinates, neglecting mismatch types. Therefore, by retaining only positions with
   SNPs from each tissue's VCF, we can compare SNPs between the two tissues and known SNPs (file 3).

I copied the file "SNPs.ClinVar.JAN24.bed.gz" to the machine, first with winscp and then with this command:

```bash
scp SNPs.ClinVar.JAN24.bed.gz maayanbah@10.0.32.173:/home/maayanbah/hw_5/SNPs.ClinVar.JAN24.bed.gz
````

```bash
cd /home/maayanbah/hw_5/clean_samples
bcftools view --types snps -o tumor_snps.vcf "tumor_samples.vcf"
bcftools view --types snps -o normal_snps.vcf "normal_samples.vcf"

bedtools intersect -v -a tumor_snps.vcf -b "/home/maayanbah/hw_5/SNPs.ClinVar.JAN24.bed.gz" normal_snps.vcf > novel_and_somatic_variants.vcf
```

The file is in the path "/home/maayanbah/hw_5/clean_samples/novel_and_somatic_variants.vcf"

## Step 6

Now we'll use Salmon to quantify per-sample expression levels.

```bash
# activate conda
conda activate compgen

# download transcriptome & index it with salmon
mkdir -p GenomicPractice/genomes/human/Transcripome/
cd GenomicPractice/genomes/human/Transcripome/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
salmon index -t refMrna.fa.gz -i hg38_refGene_index

# quantify expression with salmon
cd ~/hw_5/clean_samples
mkdir salmon_results

accessions=(
   "ERR031028"
   "ERR031030"
   "ERR031032"
   "ERR299297"
   "ERR299298"
   "ERR031038"
   "ERR031040"
   "ERR031042"
   "ERR031044"
   "ERR031018"
   "ERR299295"
   "ERR031022"
   "ERR031024"
   "ERR031027"
   "ERR031029"
   "ERR031031"
   "ERR031033"
   "ERR031035"
   "ERR299299"
   "ERR031039"
   "ERR031041"
   "ERR031043"
   "ERR031017"
   "ERR031019"
   "ERR299296"
   "ERR031023"
)

for acc in ${accessions[@]}; do echo $acc; 
 echo salmon quant -i ~/GenomicPractice/genomes/human/Transcripome/hg38_refGene_index -l A -1 ${acc}_1.fastp.fastq.gz -2 \
 ${acc}_2.fastp.fastq.gz -p 8 --validateMappings  -o salmon_results/${acc};
 salmon quant -i ~/GenomicPractice/genomes/human/Transcripome/hg38_refGene_index -l A -1 ${acc}_1.fastp.fastq.gz -2 \
 ${acc}_2.fastp.fastq.gz -p 8 --validateMappings  -o salmon_results/${acc}

done

```

Now I uploaded:

1. transcript_to_gene.py
2. TranscriptToGene.hg38.tsv

To the class machine using winscp and:

```bash
scp TranscriptToGene.hg38.tsv maayanbah@10.0.32.173:/home/maayanbah/hw_5/TranscriptToGene.hg38.tsv
scp transcript_to_gene.py maayanbah@10.0.32.173:/home/maayanbah/hw_5/transcript_to_gene.py
```

Then, we'll use the script transcript_to_gene.py to obtain one file containing the per-sample per-gene reads counts and
TPM levels.

```bash
# combine per-sample salmon's results into one file
python ~/hw_5/transcript_to_gene.py ~/hw_5/TranscriptToGene.hg38.tsv ~/hw_5/clean_samples/salmon_results \
~/hw_5/clean_samples/salmon_results/out_file_merged.tsv
```

## Step 7

We'll use PyDESeq2 to analyze differential expression between the tissues, using cutoffs of
|𝑙𝑜𝑔2(𝑓𝑜𝑙𝑑𝑐ℎ𝑎𝑛𝑔𝑒(𝑥))|≥2 and −𝑙𝑜𝑔10(𝑝𝑎𝑑𝑗(𝑥))≥−𝑙𝑜𝑔10(5∗10−5), where 𝑥 is a gene.
we will create two gene lists and save them to the disk:

1. Significantly changed genes (according to cutoffs).
2. Significantly positively changed genes, that it, genes whose expression levels increased in tumor tissues
   significantly.

The python script COVID19.DiffExp.PyDESeq2.py from the practice with minor modifications:

```python
#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np
import dash_bio
import plotly.express as px


def load_data(counts_file):
    counts_df = pd.read_table(counts_file).pivot(
        columns="Gene",
        index="Sample",
        values="NumReads"
    ).sort_index().astype(int)

    tumor_samples = (
        "ERR031028 ERR031030 ERR031032 ERR299297 ERR299298 ERR031038 ERR031040 ERR031042 ERR031044 ERR031018 ERR299295"
        " ERR031022 ERR031024"
    )
    normal_samples = (
        "ERR031027 ERR031029 ERR031031 ERR031033 ERR031035 ERR299299 ERR031039 ERR031041 ERR031043 ERR031017 ERR031019"
        " ERR299296 ERR031023"
    )

    metadata = pd.DataFrame(
        {
            "prostate_cancer": ["Noraml" for _ in range(1, 14)] + ["Tumor" for _ in range(1, 14)]},
        index=f"{normal_samples} {tumor_samples}".split(),
    )

    return counts_df, metadata


def preprocess_data(counts_df, metadata):
    """
    Preprocess data by removing samples with missing annotations and genes with low expression.
    """
    samples_to_keep = ~metadata["prostate_cancer"].isna()
    counts_df = counts_df.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]

    # Next, we filter out genes that have less than 10 read counts in total.
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]

    # Sort both counts_df and metadata based on sample indices
    counts_df = counts_df.sort_index()
    metadata = metadata.sort_index()

    return counts_df, metadata


def create_show_and_save_fig(significantly_changed_genes_df,
                             clustergram_image_name,
                             clustergram_title,
                             counts_file,
                             metadata,
                             out_dir):
    """Generate, display, and save a clustergram figure."""
    significant_tpm_df = (
        pd.read_table(counts_file)
        .pivot(
            columns="Gene",
            index="Sample",
            values="TPM"
        )
        .sort_index()
        .astype(int)
        .sort_index()
        .T
        .rename_axis(None, axis=1)
        .reset_index()
        .pipe(lambda df: df.loc[df["Gene"].isin(significantly_changed_genes_df.values)])
        .set_index("Gene")
        .T
        .reset_index()
        .rename_axis(None, axis=1)
        .merge(metadata.reset_index(), how="left")
        .rename(columns={"index": "Sample", "prostate_cancer": "Prostate Cancer"})
        .set_index(["Sample", "Prostate Cancer"])
        .T
    )

    fig = dash_bio.Clustergram(
        data=significant_tpm_df,
        column_labels=list(significant_tpm_df.columns.values),
        row_labels=list(significant_tpm_df.index),
        cluster="row",
        column_colors=["green" for _ in range(1, 14)] + ["red" for _ in range(1, 14)],
        height=400,
        width=700
    )
    fig.update_layout(title=clustergram_title)
    fig.write_image(Path(out_dir, clustergram_image_name), width=400, height=700)
    fig.show()


def generate_volcano_plot(stat_res_df, out_dir):
    """Generate and display a volcano plot."""
    stat_res_df["-log10(padj)"] = -np.log10(stat_res_df["padj"])

    fig = px.scatter(
        stat_res_df,
        x="log2FoldChange",
        y="-log10(padj)",
        title="A volcano plot demonstrating the significance of change in expression due to prostate cancer",
        hover_data=["Gene"]
    )
    fig.add_shape(
        type="line",
        x0=-1,
        x1=-1,
        y0=0,
        y1=stat_res_df["-log10(padj)"].max(),
        line=dict(
            color="red",
            width=2,
            dash="dash"
        ),
    )
    fig.add_shape(
        type="line",
        x0=1,
        x1=1,
        y0=0,
        y1=stat_res_df["-log10(padj)"].max(),
        line=dict(
            color="red",
            width=2,
            dash="dash"
        ),
    )
    fig.add_shape(
        type="line",
        x0=stat_res_df["log2FoldChange"].min(),
        x1=stat_res_df["log2FoldChange"].max(),
        y0=-np.log10(0.05),
        y1=-np.log10(0.05),
        line=dict(
            color="black",
            width=2,
            dash="dash"
        ),
    )
    fig.update_layout(template="plotly_white", width=600, height=400)
    fig.write_image(Path(out_dir, "VolacnoPlot.png"), width=600, height=400)
    fig.show()

def main():
    out_dir = Path("~/hw_5/clean_samples/salmon_results").expanduser()
    if not out_dir.is_dir():
        raise FileNotFoundError(f"{out_dir} is not a valid directory.")

    counts_file = Path("~/hw_5/clean_samples/salmon_results/out_file_merged.tsv").expanduser()
    counts_df, metadata = load_data(counts_file)

    counts_df, metadata = preprocess_data(counts_df, metadata)

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="prostate_cancer",
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()
    print(dds)

    """
    Now that dispersions and LFCs were fitted, we may proceed with statistical tests to compute p-values and adjusted
    p-values for differential expression. This is the role of the DeseqStats class. It has a unique mandatory
    argument, dds, which should be a fitted DeseqDataSet object.
    """

    stat_res = DeseqStats(dds, inference=inference)

    """
    The `summary` function runs the statistical analysis (multiple testing adjustment included) and returns a summary
    DataFrame.
    """
    stat_res.summary()
    dds.varm["LFC"]

    stat_res_df = stat_res.results_df.reset_index()
    generate_volcano_plot(stat_res_df, out_dir)

    # Filter genes based on conditions
    significantly_positively_changed_genes = stat_res_df.loc[
        (stat_res_df['log2FoldChange'] >= 2) &
        (stat_res_df["-log10(padj)"] >= -np.log10(5 * 10 ** -5)),
        "Gene"
    ]

    significantly_changed_genes = stat_res_df.loc[
        (stat_res_df["log2FoldChange"].abs() >= 2) &
        (stat_res_df["-log10(padj)"] >= -np.log10(5 * 10 ** -5)),
        "Gene"
    ]

    create_show_and_save_fig(
        significantly_changed_genes,
        "Clustergram_significantly_changed_genes.png",
        "Clustergram (clustered heatmap) of significantly changed genes",
        counts_file,
        metadata,
        out_dir
    )

    significantly_changed_genes.to_csv("~/hw_5/significantly_changed_genes.csv", index=False)
    create_show_and_save_fig(
        significantly_positively_changed_genes,
        "Clustergram_positively_significantly_changed_genes.png",
        "Clustergram (clustered heatmap) of positively significantly changed genes",
        counts_file,
        metadata,
        out_dir
    )
    # Save significantly changed genes and significantly positively changed genes to disk
    significantly_positively_changed_genes.to_csv("~/hw_5/significantly_positively_changed_genes.csv", index=False)


if __name__ == "__main__":
    main()

```

```bash
pip install pydeseq2
# Run the code
python COVID19.DiffExp.PyDESeq2.py
# Download files to my computer using scp (snd after that using WinSCP)
scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/clean_samples/salmon_results/VolacnoPlot.png VolacnoPlot.png
scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/clean_samples/salmon_results/Clustergram_positively_significantly_changed_genes.png Clustergram_positively_significantly_changed_genes.png
scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/clean_samples/salmon_results/Clustergram_significantly_changed_genes.png Clustergram_significantly_changed_genes.png
scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/significantly_positively_changed_genes.csv significantly_positively_changed_genes.csv
scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/significantly_changed_genes.csv significantly_changed_genes.csv
```

## Step 8

The volcano plot from step 7:
![image](https://github.com/MaayanBah/prostate-cancer-rna-analysis/assets/84293984/e69efa9e-4ad4-4809-96c0-36cc0e6ce679)

The clustergram of the significantly changed genes and the positively significantly changed genes:
![image](https://github.com/MaayanBah/prostate-cancer-rna-analysis/assets/84293984/701fc400-2280-4e9b-9ac9-0c0082449145)

![image](https://github.com/MaayanBah/prostate-cancer-rna-analysis/assets/84293984/67788fcb-e7a1-4f2a-9c97-201a425fcbd2)

GO analysis on the significantly changed genes:
![image](https://github.com/MaayanBah/prostate-cancer-rna-analysis/assets/84293984/6e7e6bb3-95a2-4c6d-9ca4-bf6e238320a5)

## Step 9

Finally, we want to combine our insights from novel somatic variants and differentially expressed genes in prostate
cancer. Using only positively significantly changed genes (from step 7), we will count the number of novel, somatic
variants (from step 5) per each such gene.

We will supply you with a list of the 20 genes with the most such SNPs inside, sorted by descending order of the number
of SNPs.

We'll do this as follows:

1. We'll link the gene symbols of the positively significantly changed genes to their gene ID equivalents (file 4)
   and use that to link these genes to their coordinates in the genes bed file (file 2).

   ```python
   import csv
   import pandas as pd
   
   
   def create_gene_name_to_id_dict() -> dict[str, list[str]]:
       name_to_id: dict[str, list[str]] = {}
       with open("GeneIDGeneSymbol.txt") as csvfile:
           reader = csv.reader(csvfile)
           for row in reader:
               name: str = row[1]
               gene_id: str = row[0]
               if name in name_to_id:
                   name_to_id[name].append(gene_id)
               else:
                   name_to_id[name] = [gene_id]
       return name_to_id
   
   def main():
       name_to_id: dict[str, list[str]] = create_gene_name_to_id_dict()
   
       # Read significantly_positively_changed_genes.csv and change the symbols
       positively_changed_genes_df: pd.DataFrame = pd.read_csv(
           "step_7/significantly_positively_changed_genes.csv", sep=","
       )
   
       positively_changed_genes_ids: list[str] = []
       for gene in positively_changed_genes_df["Gene"]:
           if gene not in name_to_id:
               positively_changed_genes_ids.append(gene)
           else:
               for gene_id in name_to_id[gene]:
                   positively_changed_genes_ids.append(gene_id)
   
       # Read Genes.RefSeq.JAN24.bed.gz into a DataFrame
       bed_columns: list[str] = [
           "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
           "blockSizes", "blockStarts"
       ]
       gene_bed_df: pd.DataFrame = pd.read_csv(
           "Genes.RefSeq.JAN24.bed.gz", sep="\t", header=None, names=bed_columns
       )
   
       # Filter gene_bed_df to keep only rows where the gene ID is in positively_changed_genes_ids
       matching_genes_df: pd.DataFrame = gene_bed_df[
           gene_bed_df["name"].isin([gene_id for gene_id in positively_changed_genes_ids])
       ]
   
    matching_genes_df.to_csv("positively_changed_genes.bed", sep="\t", header=False, index=False)
   
   
   if __name__ == "__main__":
       main()
   
   ```

2. We will upload the file to the machine:

   ```bash
   scp positively_changed_genes.bed maayanbah@10.0.32.173:/home/maayanbah/hw_5/positively_changed_genes.bed
   ```
   And then we'll sort the gene coordinates from (1) and merge them by running
   bedtools merge -i <path_to_file> -s -c 4,5,6 -o distinct, sum, distinct.

   ```bash
   bedtools sort -i positively_changed_genes.bed > positively_changed_genes_sorted.bed
   
   # -s: This flag signifies that the input intervals are sorted, enhancing the merging process's efficiency.
   # -c 4,5,6: This parameter designates the columns for merging: name, score, and strand of each interval.
   # These columns will be used to determine whether intervals can be merged.
   # -o distinct,sum,distinct: This parameter outlines the operations to execute on each column during merging.
   # For the name and strand columns (columns 4 and 6), only distinct values are preserved.
   # For the score column (column 5), values are summed.
   bedtools merge -i positively_changed_genes_sorted.bed -s -c 4,5,6 -o distinct,sum,distinct > positively_changed_genes_sorted_merged.bed
   ```
   The second command combines overlapping gene regions. It combines, counts, and summarizes specific features within
   these regions, such as the distinct values in columns 4, 5, and 6, and saves the merged data into a new file
   called "positively_changed_genes_sorted_merged.bed".

3. Now we will intersect the merged gene coordinates from (2) with the novel somatic variants VCF file and count the
   number of SNPs per interval.

   ```bash
   cat /home/maayanbah/hw_5/clean_samples/novel_and_somatic_variants.vcf | grep -v \# | awk '{print $1 "\t" $2 - 1 "\t" $2}' > /home/maayanbah/hw_5/clean_samples/novel_and_somatic_variants.bed

   bedtools intersect -a positively_changed_genes_sorted_merged.bed \
   -b /home/maayanbah/hw_5/clean_samples/novel_and_somatic_variants.bed \
   -wa -c > intersected_positively_changed_genes_with_snps.txt
   ```

   The -c option in bedtools intersect is used to report the number of occurrences of each entry in the first input file
   (-a) that overlaps with entries in the second input file (-b). his count represents the number of SNPs per interval
   in the merged gene coordinates.

4. The output of (3) is the sum of SNPs per overlapping gene IDs, I copied it to my computer:
   ```bash
   scp maayanbah@10.0.32.173:/home/maayanbah/hw_5/intersected_positively_changed_genes_with_snps.txt \
   intersected_positively_changed_genes_with_snps.txt
   ```
   In each line, we will remap the gene IDs from their gene symbol. Ideally, whenever gene IDs' coordinates overlap,
   they are all isoforms of the same gene, marked by a gene symbol. However, that's not always the case,
   so we need to remove such lines.

   ```python
   import csv
   import pandas as pd
   from pprint import pprint
   
   def create_id_to_gene_name_dict() -> dict[str, str]:
       id_to_name: dict[str, str] = {}
       with open("GeneIDGeneSymbol.txt") as csvfile:
           reader = csv.reader(csvfile)
           for row in reader:
               name: str = row[1]
               gene_id: str = row[0]
               id_to_name[gene_id] = name
       return id_to_name
   
   def main():
       id_to_name: dict[str, str] = create_id_to_gene_name_dict()
   
       intersected_positively_changed_genes_with_snps = "step_9/intersected_positively_changed_genes_with_snps.txt"
       intersected_columns: list[str] = [
           "chrom", "start", "end", "name", "score", "strand", "snps_count"
       ]
       intersected_positively_changed_genes_with_snps_df: pd.DataFrame = pd.read_csv(
           "step_9/intersected_positively_changed_genes_with_snps.txt", sep="\t", header=None, names=intersected_columns
       )
   
       intersected_positively_changed_genes_with_snps_df = intersected_positively_changed_genes_with_snps_df[
           intersected_positively_changed_genes_with_snps_df["name"].apply(lambda x: len({id_to_name[gene_id] for gene_id in x.split(",")})) == 1]
       intersected_positively_changed_genes_with_snps_df["name"] = (
           intersected_positively_changed_genes_with_snps_df["name"].apply(lambda x: id_to_name[x.split(",")[0]])
       )
   
       best_intersected_positively_changed_genes_with_snps_df = (
           intersected_positively_changed_genes_with_snps_df.sort_values(by='snps_count', ascending=False).head(20)
       )
       best_intersected_positively_changed_genes_with_snps_df.to_csv(
           r"step_9\best_20_lines.csv", sep="\t", index=False
       )
   
   
   if __name__ == "__main__":
       main()
   
   ```


The gene with the most variants is SCHLAP1.
I searched for this gene in the context of prostate cancer and found the following article:

###### Kidd SG, Carm KT, Bogaard M, Olsen LG, Bakken AC, Løvf M, Lothe RA, Axcrona K, Axcrona U, Skotheim RI. High expression of SCHLAP1 in primary prostate cancer is an independent predictor of biochemical recurrence, despite substantial heterogeneity. Neoplasia. 2021 Jun

The article suggest that elevated levels of SCHLAP1 indicate a higher chance of prostate cancer returning after
treatment, as revealed by statistical analyses. Additionally, high SCHLAP1 levels are linked to worse characteristics
of the cancer, like its grade and stage, as well as specific tissue features. Overall, finding high SCHLAP1 levels in
a cancer sample suggests a poorer outlook. Notably, this is the first time SCHLAP1 has been linked to a specific
aggressive tissue feature. Since SCHLAP1 levels vary widely, it's important to test multiple samples to accurately
assess a patient's status.

The connection between increased expression levels of SCHLAP1 and the creation of somatic mutations in SCHLAP1 itself
may be:

1. Gene regulation: SCHLAP1 can regulate gene expression through various mechanisms,
   including chromatin remodeling (1), epigenetic modifications(2) and transcriptional regulation. Increased expression
   of SCHLAP1 might disrupt these mechanisms, leading to abnormal gene expression patterns, including mutations
   within the SCHLAP1 gene itself.
2. Genomic instability: The (3) article suggests a potential link between increased expression of SCHLAP1 and
   genomic instability, which could contribute to the aggressiveness of prostate cancer subtypes characterized by
   IDC and CA. Essentially, high levels of SCHLAP1 expression may contribute to genomic instability, and Genomic
   instability increases the likelihood of mutations occurring throughout the genome, including within the SCHLAP1
   gene locus.
3. Cellular proliferation and survival: SCHLAP1 overexpression is associated with enhanced cellular proliferation
   and survival (4), which can increase the likelihood of mutations occurring due to errors during DNA
   replication or exposure to genotoxic stressors.

###### (1) Raab JR, Smith KN, Spear CC, Manner CJ, Calabrese JM, Magnuson T. SWI/SNF remains localized to chromatin in the presence of SCHLAP1. Nat Genet. 2019 Jan

###### (2) Anirban Sahu; John R. Prensner; Qi Cao; Arul M. Chinnaiyan. SChLAP1 mediated epigenetic modifications in prostate cancer. Cancer Res (2015)

###### (3) Chua MLK, Lo W, etc. A Prostate Cancer "Nimbosus": Genomic Instability and SChLAP1 Dysregulation Underpin Aggression of Intraductal and Cribriform Subpathologies. Eur Urol. 2017 Nov

###### (4) Huang, K., Tang, Y. SChLAP1 promotes prostate cancer development through interacting with EZH2 to mediate promoter methylation modification of multiple miRNAs of chromosome 5 with a DNMT3a-feedback loop. Cell Death Dis 12, 188 (2021).
