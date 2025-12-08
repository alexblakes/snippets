# Matplotlib snippets

Flatten an array of Axes

```python
# flatten() is an np.ndarray method
axs = axs.flatten() # Column-wise
axs = axs.flatten("F") # Row-wise
```

Subset a colormap
```python
import numpy as np
from matplotlib import cm
import matplotlib.colors as mcolors

def subset_colormap(lo=0.4, hi=1, colormap=cm.Reds, name="new_colormap"):
    color_space = np.linspace(lo, hi, 10)
    colors = colormap(color_space)
    return mcolors.LinearSegmentedColormap.from_list(name, colors)
```
