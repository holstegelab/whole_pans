# Build-aware annotation of only the technically passing frozen panel.
SGP_ANN = SGP.get("annotation", {})
SGP_ANN_ENV = str(SGP_ROOT / "workflow/envs/svgp_annotation.yaml")
SGP_ANN_DIR = SGP_OUT + "/annotation"
SGP_RM = SGP_ANN.get("repeatmasker", {})


rule annotation:
    input:
        SGP_ANN_DIR + "/allele_gene_context.tsv",
        SGP_ANN_DIR + "/feature_evidence.tsv",
        SGP_ANN_DIR + "/descriptive_counts.tsv"


rule svgp_annotation_prepare:
    input:
        panel=SGP_OUT + "/panel/panel.tsv",
        fasta=SGP_OUT + "/catalogue/alleles.fa",
        lock=SGP["measurement"]["measurement_lock"],
        features=[SGP_ANN["features"]] if SGP_ANN.get("features") else [],
        review=[SGP_ANN["coordinate_review"]] if SGP_ANN.get("coordinate_review") else []
    output:
        loci=SGP_ANN_DIR + "/passing_loci.tsv",
        metadata=SGP_ANN_DIR + "/features.tsv",
        features=SGP_ANN_DIR + "/features.raw.bed",
        genes=SGP_ANN_DIR + "/gene.raw.bed",
        tss=SGP_ANN_DIR + "/tss.raw.bed",
        queries=SGP_ANN_DIR + "/queries.raw.bed",
        inserted=SGP_ANN_DIR + "/inserted.fa"
    params:
        task="prepare",
        main_graph=SGP["measurement"]["main_graph"],
        settings=SGP_ANN
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    log: SGP_OUT + "/logs/annotation_prepare.log"
    conda: SGP_ANN_ENV
    script: "../scripts/svgp_annotation.py"


rule svgp_annotation_sort:
    input:
        features=SGP_ANN_DIR + "/features.raw.bed",
        genes=SGP_ANN_DIR + "/gene.raw.bed",
        tss=SGP_ANN_DIR + "/tss.raw.bed",
        queries=SGP_ANN_DIR + "/queries.raw.bed"
    output:
        features=SGP_ANN_DIR + "/features.bed",
        genes=SGP_ANN_DIR + "/gene.bed",
        tss=SGP_ANN_DIR + "/tss.bed",
        queries=SGP_ANN_DIR + "/queries.bed",
        version=SGP_ANN_DIR + "/bedtools.version.txt"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=30
    log: SGP_OUT + "/logs/annotation_sort.log"
    conda: SGP_ANN_ENV
    shell:
        """
        {{
            bedtools --version > {output.version:q}
            LC_ALL=C sort -k1,1 -k2,2n {input.features:q} > {output.features:q}
            LC_ALL=C sort -k1,1 -k2,2n {input.genes:q} > {output.genes:q}
            LC_ALL=C sort -k1,1 -k2,2n {input.tss:q} > {output.tss:q}
            LC_ALL=C sort -k1,1 -k2,2n {input.queries:q} > {output.queries:q}
        }} > {log:q} 2>&1
        """


rule svgp_annotation_intersect:
    input:
        queries=SGP_ANN_DIR + "/queries.bed",
        features=SGP_ANN_DIR + "/features.bed"
    output: SGP_ANN_DIR + "/overlaps.bed"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=30
    log: SGP_OUT + "/logs/annotation_intersect.log"
    conda: SGP_ANN_ENV
    shell:
        """
        if [ -s {input.queries:q} ] && [ -s {input.features:q} ]; then
            bedtools intersect -sorted -wa -wb -a {input.queries:q} -b {input.features:q} > {output:q} 2> {log:q}
        else
            : > {output:q}
            : > {log:q}
        fi
        """


rule svgp_annotation_nearest:
    input:
        queries=SGP_ANN_DIR + "/queries.bed",
        features=SGP_ANN_DIR + "/{kind}.bed"
    output: SGP_ANN_DIR + "/nearest_{kind}.bed"
    wildcard_constraints:
        kind="gene|tss"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=30
    log: SGP_OUT + "/logs/annotation_nearest_{kind}.log"
    conda: SGP_ANN_ENV
    shell:
        """
        if [ -s {input.queries:q} ] && [ -s {input.features:q} ]; then
            bedtools closest -a {input.queries:q} -b {input.features:q} -t all -d > {output:q} 2> {log:q}
        else
            : > {output:q}
            : > {log:q}
        fi
        """


rule svgp_annotation_tss_window:
    input:
        queries=SGP_ANN_DIR + "/queries.bed",
        tss=SGP_ANN_DIR + "/tss.bed"
    output: SGP_ANN_DIR + "/tss_within_100kb.bed"
    threads: 1
    resources:
        mem_mb=4000,
        runtime_min=30
    log: SGP_OUT + "/logs/annotation_tss_window.log"
    conda: SGP_ANN_ENV
    shell:
        """
        if [ -s {input.queries:q} ] && [ -s {input.tss:q} ]; then
            bedtools window -w 100000 -a {input.queries:q} -b {input.tss:q} > {output:q} 2> {log:q}
        else
            : > {output:q}
            : > {log:q}
        fi
        """


if SGP_RM.get("enabled", False):
    if not SGP_RM.get("library_release"):
        raise ValueError("Record annotation.repeatmasker.library_release before enabling RepeatMasker")
    rule svgp_annotation_repeatmasker:
        input:
            sequence=SGP_ANN_DIR + "/inserted.fa",
            library=SGP_RM["library"]
        output:
            repeats=SGP_ANN_DIR + "/repeatmasker/inserted.fa.out",
            version=SGP_ANN_DIR + "/repeatmasker/version.txt",
            library_checksum=SGP_ANN_DIR + "/repeatmasker/library.sha256"
        params:
            executable=SGP_RM.get("executable", "RepeatMasker"),
            out=SGP_ANN_DIR + "/repeatmasker"
        # RMBlast uses four cores per -pa batch, including -pa 1.
        threads: 4 * int(SGP_RM.get("parallel_batches", 1))
        resources:
            mem_mb=SGP_RM.get("mem_mb", 16000),
            runtime_min=SGP_RM.get("runtime_min", 240)
        log: SGP_OUT + "/logs/annotation_repeatmasker.log"
        benchmark: SGP_OUT + "/benchmarks/annotation_repeatmasker.tsv"
        conda: SGP_RM["conda_env"]
        shell:
            """
            {{
                test -s {input.sequence:q}
                test {threads} -ge 4
                {params.executable:q} -version > {output.version:q}
                sha256sum {input.library:q} > {output.library_checksum:q}
                {params.executable:q} -engine rmblast -pa $(( {threads} / 4 )) -lib {input.library:q} -dir {params.out:q} {input.sequence:q}
            }} > {log:q} 2>&1
            """


rule svgp_annotation_summary:
    input:
        loci=SGP_ANN_DIR + "/passing_loci.tsv",
        metadata=SGP_ANN_DIR + "/features.tsv",
        overlaps=SGP_ANN_DIR + "/overlaps.bed",
        nearby=SGP_ANN_DIR + "/tss_within_100kb.bed",
        nearest_gene=SGP_ANN_DIR + "/nearest_gene.bed",
        nearest_tss=SGP_ANN_DIR + "/nearest_tss.bed",
        links=[SGP_ANN["regulatory_links"]] if SGP_ANN.get("regulatory_links") else [],
        repeatmasker=[SGP_ANN_DIR + "/repeatmasker/inserted.fa.out"] if SGP_RM.get("enabled", False) else []
    output:
        context=SGP_ANN_DIR + "/allele_gene_context.tsv",
        evidence=SGP_ANN_DIR + "/feature_evidence.tsv",
        counts=SGP_ANN_DIR + "/descriptive_counts.tsv"
    params:
        task="summarize",
        settings=SGP_ANN
    threads: 1
    resources:
        mem_mb=8000,
        runtime_min=60
    log: SGP_OUT + "/logs/annotation_summary.log"
    conda: SGP_ANN_ENV
    script: "../scripts/svgp_annotation.py"
