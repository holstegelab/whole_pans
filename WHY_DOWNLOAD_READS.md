# Why Download Reads Instead of Using Assemblies Alone?

## Main reason

Reads provide **independent evidence** for structural variants discovered from assemblies. An assembly is a reconstruction of the reads, so assembly errors can resemble novel variants. Using the same assembly for both discovery and validation does not establish that the allele is real.

Assembly-derived false positives can result from:

- collapsed or expanded repeats;
- chimeric contig joins;
- haplotype switches;
- polishing errors;
- fragmented contigs;
- missing heterozygous alleles;
- inaccurate sequence or breakpoints in repetitive regions.

Read alignments can show whether multiple independent molecules support or span a candidate breakpoint. They also help distinguish a true reference genotype from an assembly omission or a low-coverage no-call.

## Role of each data type

- **Assemblies:** discover candidate alleles and provide candidate sequences and haplotype context.
- **Reads:** validate candidates, refine breakpoints, genotype samples, measure allele balance, and distinguish absence from missing information.

Assemblies are therefore sufficient for a **provisional discovery catalog**, but reads are needed before uncertain alleles are accepted into a high-confidence catalog or a new graph release. Read-level genotyping is also necessary for reliable allele-frequency and population analyses.

## Relevance to the current workflow

The current graph screen shows that assemblies with poor contiguity can produce inflated candidate counts despite acceptable graph callability. Sample `5352834` is a notable example and should be treated as assembly-driven until its candidates receive independent support.

Agreement between graph-based discovery, a linear-reference assembly caller, and multiple unrelated samples increases confidence. However, singleton, complex, repetitive, or poor-assembly-only candidates remain especially dependent on read validation.

## Practical download priority

It is not necessary to download reads for every sample immediately. Prioritize:

1. carriers of high-value candidate alleles;
2. fragmented or failed assemblies;
3. singleton and complex structural variants;
4. candidates supported by multiple discovery methods but not yet by reads;
5. graph-member samples used as negative or positive controls.

For complex structural variants and repeats, long accurate reads such as PacBio HiFi are preferable. Short reads can support genotyping of represented alleles, but they may not resolve complex breakpoints or long repetitive insertions.

## Decision rule

> Use assemblies to discover candidate variants, but require independent read support—or exceptionally strong orthogonal and multi-sample evidence—before treating an uncertain allele as validated or adding it to graph v2.
