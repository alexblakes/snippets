# Python snippets

## Dedent multiline strings
Use [`inspect.cleandoc()`](https://docs.python.org/3/library/inspect.html#inspect.cleandoc)

```python
vcf_header = inspect.cleandoc(
  """
  ##fileformat=VCFv4.2
  ##contig=<ID=chr11,length=135086622>
  #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
  """
)
```

See [this discussion](https://discuss.python.org/t/indented-multi-line-string-literals/9846/2) at python.org.

See also [`textwrap.dedent()`](https://docs.python.org/3/library/textwrap.html#textwrap.dedent)

Or the utility function described [here](https://discuss.python.org/t/indented-multi-line-string-literals/9846/2).
