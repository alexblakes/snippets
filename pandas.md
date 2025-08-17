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

## Methods chains to assign with functions that return multiple values
I am often running statistical tests row-wise. These tests are defined in functions which return multiple values (e.g. `odds_ratio`, `ci_lo`, `ci_hi`, `p`). The method chaining syntax for this is hard to remember. Here's a generic function which can be used in place.
```python
def assign_with_per_row_fn(df, fn, new_cols, *args, **kwargs):
    return df.assign(
        **pd.DataFrame(
        df.apply(fn, axis=1, args=tuple(*args), **kwargs).to_list(),
        columns=new_cols,
        index=df.index
        )
    )
```
