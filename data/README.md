# AVA,Dx data sources

## Sources available online - downloaded automatically:
```
.
├── annovar
│   └── humandb
│       ├── annovar_downdb_gnomad_exome.log
│       ├── annovar_downdb_gnomad_genome.log
│       ├── annovar_downdb_refGene.log
│       ├── <assembly_version>_gnomad_exome.txt.gz
│       ├── <assembly_version>_gnomad_exome.txt.idx
│       ├── <assembly_version>_gnomad_genome.txt.gz
│       ├── <assembly_version>_gnomad_genome.txt.idx
│       ├── <assembly_version>_refGene.txt
│       ├── <assembly_version>_refGeneMrna.fa
│       └── <assembly_version>_refGeneVersion.txt
├── avadx
│   ├── CPDB_pathways_genesymbol.tab
│   ├── varidb.db
│   ├── varidb.log
│   └── varidb.md5
├── ethseq
│   └── models
│       └── Exonic.All.Model.gds
└── refseq
    ├── <assembly_version>_feature_table.txt
    └── <assembly_version>_protein.faa
```

## Sources created based on downloaded data:
```
.
└── avadx
    ├── Transcript-ProtLength.csv
    ├── Transcript-ProtLength_cleaned.csv
    ├── <assembly_version>_gnomad_exome_allAFabove0.txt.gz
    ├── <assembly_version>_gnomad_exome_allAFabove0.txt.gz.tbi
    ├── <assembly_version>_gnomad_genome_allAFabove0.txt.gz
    ├── <assembly_version>_gnomad_genome_allAFabove0.txt.gz.tbi
    ├── prot_seqs.fa
    └── refseq_mapping.csv
```
