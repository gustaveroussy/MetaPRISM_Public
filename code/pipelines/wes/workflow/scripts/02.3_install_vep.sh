#!/bin/bash

while getopts ":a:s:g:c:p:f:r:" opt; do
    case $opt in
        a) assembly="$OPTARG"
            ;;
        s) species="$OPTARG"
            ;;
        g) plugins="$OPTARG"
            ;;
        c) cache="$OPTARG"
            ;;
        p) path="$OPTARG"
            ;;
        f) fasta="$OPTARG"
            ;;
        r) release="$OPTARG"
            ;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac

    case $OPTARG in
        -*) echo "Option $opt needs a valid argument"
            exit 1
            ;;
    esac
done

# Load modules
module load gcc

# message
printf "working directory: %s\n" $PWD

# Get full paths before changing directory
cache_full_path=$PWD/${cache}

# Get the code
if [[ -d ${path} ]]
then
    rm -r ${path}
fi

curl -L -o ${path}.zip https://github.com/Ensembl/ensembl-vep/archive/release/${release}.0.zip
unzip -qq ${path}.zip
mv ensembl-vep-release-${release}.0 ${path}
rm ${path}.zip
cd ${path}

# Install perl libraries
perl_libs=$PWD/perl_libs
mkdir -p ${perl_libs}

printf "INSTALLING Bio::SeqFeature::Lite ...\n"
cpanm -l ${perl_libs} Bio::SeqFeature::Lite

printf "INSTALLING Bio::DB::HTS::Tabix ...\n"
cpanm -l ${perl_libs} Bio::DB::HTS::Tabix

export PERL5LIB=${perl_libs}/lib/perl5:$PERL5LIB

perl INSTALL.pl --AUTO afc \
    --NO_UPDATE \
    --NO_HTSLIB \
    --ASSEMBLY ${assembly} \
    --SPECIES ${species} \
    --PLUGINS ${plugins} \
    --CACHEDIR ${cache_full_path}

# Uncompressing the fasta
gunzip ${fasta}
