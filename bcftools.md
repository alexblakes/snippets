# bcftools snippets and gotchas

## Piping with `bcftools annotate`

VCF annotation files seem not to be supported when input is piped into `bcftools annotate`.

Annotation files in other formats, e.g. TSV, do seem to be supported.

I'm not clear whether this is a bug.

## `bcftools merge` creates multi-allelic sites

When VCFs are merged with `bcftools merge`, alleles at overlapping positions are merged into multiallelic sites.

Use the `--merge none` flag to prevent this behaviour.

## `bcftools +fill-tags`

The fill-tags plugin will throw an error if the `--sample-file` (`-S`) contains duplicates.

## `samtools --targets-file` and `--regions-file`

The `bcftools -T` / `--targets-file` option takes a bed file and limits your operation to regions within the bed file. However, it iterates through the entire VCF checking for overlap with this bed file. For large VCFs this can take a long time.

The `bcftools -R` / `--regions-file` option does the same thing, but uses the VCF's index to hop directly to the region of interest. This saves you iterating through the whole VCF. However, there is a small overhead to retrieving the index, so if you have a huge number of regions in your bed file, it may be quicker to use -T instead.

(`-T` and `-R` can also be used jointly)

In samtools it looks as though the `-L` / `--targets-file` and the `--regions-file` behave similarly to bcftools. It also has an optional region argument to be specified on the command line. I found the documentation for these was much more opaque.

A bit of trial and error meant a massive speed up using samtools view --regions-file example.bed example.bam chr1 for one particular use case.
