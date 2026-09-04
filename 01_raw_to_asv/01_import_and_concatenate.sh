#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 01.1 — Concatenate lanes per sample and import into QIIME2
# ══════════════════════════════════════════════════════════════════════════
# Each sample was sequenced across two lanes (four FASTQ files total: R1/R2
# x lane 1/lane 2). QIIME2 expects one R1 + one R2 per sample, so lanes are
# concatenated first.

mkdir -p ../02.qiime_ready

for d in */ ; do
  sample=${d%/}
  echo "Procesando $sample"

  r1_files=$(ls "$sample"/*_L*_1.fq.gz 2>/dev/null | sort)
  r2_files=$(ls "$sample"/*_L*_2.fq.gz 2>/dev/null | sort)

  n1=$(echo "$r1_files" | wc -w)
  n2=$(echo "$r2_files" | wc -w)

  if [ "$n1" -eq 0 ] || [ "$n2" -eq 0 ]; then
    echo "  No encuentro los fastq en $sample (revisa nombres o subcarpetas)"
    continue
  fi

  cat $r1_files > "../02.qiime_ready/${sample}_S1_L001_R1_001.fastq.gz"
  cat $r2_files > "../02.qiime_ready/${sample}_S1_L001_R2_001.fastq.gz"

  echo "  R1: $n1 archivos | R2: $n2 archivos"
done

cd ../02.qiime_ready
ls *_R1_001.fastq.gz | wc -l
ls *_R2_001.fastq.gz | wc -l
cd ..

# ---- Import into QIIME2 ----
conda activate qiime2-2022.11

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 02.qiime_ready \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path 02.qiime_ready/demux-paired-end.qza

qiime demux summarize \
  --i-data 02.qiime_ready/demux-paired-end.qza \
  --o-visualization 02.qiime_ready/demux-paired-end.qzv
