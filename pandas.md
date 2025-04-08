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
column_name = "my_column"
df = df.assign(**{column_name: df['A'] + df['B']})
```

## Flatten heirarchical columns
```python
def flatten_columns(df):
    df.columns = ["_".join(str(x) for x in tuple) for tuple in df.columns.to_flat_index()]
    return df

df.pipe(flatten_columns)
```
