# Matplotlib snippets

Flatten an array of Axes

```python
# flatten() is an np.ndarray method
axs = axs.flatten() # Column-wise
axs = axs.flatten("F") # Row-wise
```