
WKDIR = '/mnt/home/guojiaro/Documents/software/gits/uniqprimer/uniqprimer/test/data'
SAMPLES = 'mock'
DBDIR='/mnt/research/tg/g/glbrc_my/contam/jgi/torben_asse/db' #att: use abspath
SCRIPTDIR='/mnt/research/tg/g/glbrc_my/contam/jgi/torben_asse/script'
TOPN=2 # pick top n most abundant genome bins as target
K=20

workdir: WKDIR

localrules:

rule all:
    input:
        expand(
            [directory('{sample}/maxbin'),
            '{sample}/maxbin/maxbin.summary',
            '{sample}/mb.cm.out',
            '{sample}/mb.cm.tax',
            '{sample}/mb.mashscreen.tsv',
            '{sample}/mb.mashscreenw.tsv',
            '{sample}/maxbin/maxbin.{idx}.fasta.mashdist.tsv',
            '{sample}/mb.rmcontain.csv',
            '{sample}/maxbin/maxbin.{idx}.fasta.rmmatch.csv',
            '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer.uniq2ref.primer'],
            sample=SAMPLES, idx=['{:0>3d}'.format(i) for i range(1, TOPN+1)]
        )

rule maxbin:
    input: 
        '{sample}/{pid}_contigs.fasta.gz',
        '{sample}/{pid}_contigcoverage.tsv'
    output: 
        directory('{sample}/maxbin'),
        '{sample}/maxbin/maxbin.summary'
    conda: 'envs/root.yaml'
    threads: 4
    shell: """
        cd {wildcards.sample} \
        && mkdir -p maxbin \
        && run_MaxBin.pl -thread {threads} \
            -max_iteration 100 \
            -contig {wildcards.pid}_contigs.fasta.gz \
            -abund {wildcards.pid}_contigcoverage.tsv \
            -out maxbin/maxbin
    """
        
rule checkm:
    input:
        expand(
            '{{sample}}/maxbin/maxbin.{idx}.fasta',
            idx=['{:0>3d}'.format(i) for i in range(1, TOPN+1)]
        )
    output:
        '{sample}/mb.cm.out',
        '{sample}/mb.cm.tax'
    conda: 'envs/checkm.yaml'
    threads: 4
    shell: """
        cd {wildcards.sample}
        rm -rf mb.cm.tmpdir
        checkm lineage_wf -x fasta \
            -f mb.cm.out \
            -t {threads} \
            maxbin \
            mb.cm.tmpdir \
          || rm -rf mb.cm.tmpdir
        checkm tree_qa mb.cm.tmpdir -o 1 -f mb.cm.tax --tab_table
    """

rule mash_screen:
    input:
        '{sample}/{sample}_reads.fasta.gz'
    output:
        '{sample}/mb.mashscreen.tsv',
        '{sample}/mb.mashscreenw.tsv'
    conda: 'envs/root.yaml'
    threads: 4
    shell: """
        cd {wildcards.sample}
        ## gather
        mash screen -p {threads} {DBDIR}/refseq.genomes.k21s1000.msh {input} \
          | sort -g -r > mb.mashscreen.tsv
        mash screen -w -p {threads} {DBDIR}/refseq.genomes.k21s1000.msh {input} \
          | sort -g -r > mb.mashscreenw.tsv
    """
        
rule mash_dist:
    input:
        '{sample}/maxbin/maxbin.{idx}.fasta'
    output:
        '{sample}/maxbin/maxbin.{idx}.fasta.mashdist.tsv'
    conda: 'pipeline/envs/root.yaml'
    threads: 4
    shell: """
        # only pick top 10 most abundant to analyze
        if [ {wildcards.idx} -gt {TOPN} ]; then exit 0; fi
        cd {wildcards.sample}/maxbin
        ## search
        # -m: minimum kmer count, 1 for genome or bins, >2 for reads
        f=maxbin.{wildcards.idx}.fasta
        mash sketch -p {threads} -k 21 -m 1 $f
        # sort -k 3 --> pick 3rd col as key for sorting
        mash dist -p {threads} {DBDIR}/refseq.genomes.k21s1000.msh $f.msh \
          | sort -g -k 3 \
          > $f.mashdist.tsv
    """

# same as mash_screen with more detailed taxon info
rule refseqmasher_contain:
    input:
        '{sample}/{sample}_reads.fasta.gz'
    output:
        '{sample}/mb.rmcontain.csv'
    conda: 'envs/root.yaml'
    threads: 4
    shell: """
        cd {wildcards.sample}
        ## gather
        refseq_masher contains {input} \
        -o mb.rmcontain.csv \
        --output-type csv \
        -i 0.5 \
        -p {threads}
    """

rule refseqmasher_matches:
    input:
        '{sample}/maxbin/maxbin.{idx}.fasta'
    output:
        '{sample}/maxbin/maxbin.{idx}.fasta.rmmatch.csv'
    conda: 'envs/root.yaml'
    threads: 4
    shell: """
        # only pick top 10 most abundant to analyze
        if [ {wildcards.idx} -gt {TOPN} ]; then exit 0; fi
        cd {wildcards.sample}/maxbin
        ## search
        refseq_masher matches maxbin.{wildcards.idx}.fasta \
            -o maxbin.{wildcards.idx}.fasta.rmmatch.csv \
            --output-type csv \
            -p {threads} \
            -n 10 \ top n matches
            -m 1  # minimum kmer count, 1 for genome or bins, >2 for reads
    """

rule kmeroverlap_otherbins:
    input:
        '{sample}/maxbin/maxbin.{idx}.fasta'
    output:
        '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer'
    conda: 'env/sourmash.yaml'
    shell: """
        if [ {wildcards.idx} -gt {TOPN} ]; then exit 0; fi
        cd {wildcards.sample}/maxbin
        Seqfile_list=(maxbin.???.fasta)
        f=maxbin.{wildcards.idx}.fasta
        Other=${{Seqfile_list[@]/${{f}}}}
        python {SCRIPTDIR}/get-uniqkmer.py $f \
            --ref $Other \
            -k {K} \
            -x 4e9 -N 4 \
            --x2 1e8 --N2 4 \
          > $f.uniqkmer \
          || rm $f.uniqkmer
    """

rule kmeroverlap_combinedrefseq:
    input:
        expand(
            '{{sample}}/maxbin/maxbin.{idx}.fasta.uniqkmer',
            idx=['{:0>3d}'.format(i) for i in range(1, TOPN+1)]
        )
    output:
        '{sample}/maxbin/combined_uniqkmer.uniq2ref',
    conda: 'env/sourmash.yaml'
    shell: """
        cat {input} \
          | python {SCRIPTDIR}/get-uniqkmer-nothread.py - \
            --load {DBDIR}/refseq.20.bf \
            -k {K} --x2 1e6 --N2 4  \
          > {wildcards.sample}/maxbin/combined_uniqkmer.uniq2ref \
          || rm {wildcards.sample}/maxbin/combined_uniqkmer.uniq2ref

    """

rule kmeroverlap_persamplerefseq:
    input:
        '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer',
        '{sample}/maxbin/combined_uniqkmer.uniq2ref'
    output:
        '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer.uniq2ref',
    conda: 'env/sourmash.yaml'
    shell: """
        cd {wildcards.sample}/maxbin
        f=maxbin.{wildcards.idx}.fasta.uniqkmer
        python {SCRIPTDIR}/get-uniqkmer-nothread.py $f \
            --ref combined_uniqkmer.uniq2ref \
            -k {K} -x 1e6 -N 4 \
            --x2 1e6 --N2 4  \
          >  $f.uniq2ref \
          || rm $f.uniq2ref 
    """


rule uniqkmer_to_primer:
    input:
        '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer.uniq2ref',
        '{sample}/maxbin/maxbin.{idx}.fasta'
    output:
        '{sample}/maxbin/maxbin.{idx}.fasta.uniqkmer.uniq2ref.primer'
    conda: 'env/sourmash.yaml'
    shell: """
        cd {wildcards.sample}/maxbin
        f=maxbin.{wildcards.idx}.fasta
        echo "*** output results for $f"
        python  $Scriptdir/kmer2primer.py $Scriptdir/params.config \
          $f $f.uniqkmer.uniq2ref \
          | python $Scriptdir/filter-primer-individual.py \
            $Scriptdir/params.config - \
          | python $Scriptdir/filter-primer-pair.py $Scriptdir/params.config - \
          | tee $f.uniqkmer.uniq2ref.primer
    """
