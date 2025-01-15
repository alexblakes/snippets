# Pandas snippets

## Query on column names with spaces
Surround the column name with backticks:
```
df.query('`Genome Build` == "GRCh38"')
```
See the docs here: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html