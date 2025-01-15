# bcftools snippets and gotchas

## Piping with `bcftools annotate`

VCF annotation files seem not to be supported when input is piped into `bcftools annotate`.

Annotation files in other formats, e.g. TSV, do seem to be supported.

I'm not clear whether this is a bug.

## `bcftools merge` creates multi-allelic sites

When VCFs are merged with `bcftools merge`, alleles at overlapping positions are merged into multiallelic sites.

Use the `--merge none` flag to prevent this behaviour.