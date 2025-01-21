# Pandas snippets

## Query on column names with spaces
Surround the column name with backticks:
```
df.query('`Genome Build` == "GRCh38"')
```
See the docs here: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html

## Dynamic column names when chaining `assign()`
Use dictionary unpacking
```python
column_name = my_column
df = df.assign(**{column_name: df['A'] + df['B']})
```
