import gzip
import os

files = {
    '1m': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138515_RL2367_1m_snATAC_fragments.tsv.gz',
    '3m': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138518_RL1914_3m_snATAC_fragments.tsv.gz',
    '6m': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138521_RL2208_6m_snATAC_fragments.tsv.gz',
    '10m': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138523_RL2371_10m_snATAC_fragments.tsv.gz',
    '1y': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138526_RL2209_1y_snATAC_fragments.tsv.gz',
    '2y': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138529_RL1784_2y_snATAC_fragments.tsv.gz',
    '4y': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138532_RL2210_4y_snATAC_fragments.tsv.gz',
    'ga22': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138510_RL2366_ga22_snATAC_fragments.tsv.gz',
    'ga24': '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/GSM5138512_RL2207_ga24_snATAC_fragments.tsv.gz'
}

output_file = '/group/testa/michal.kubacki/herring/data/GSE168408_RAW/merged_output.tsv.gz'

with gzip.open(output_file, 'wt') as outfile:
    for key, file_path in files.items():
        with gzip.open(file_path, 'rt') as infile:
            for line in infile:
                outfile.write(line)
                outfile.write('\t' + key + '\n')  # Append the key to each linex``