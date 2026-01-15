## Conda in GEL

### BlockingIOError [Errno: 11]. 
Conda throws an unexpected error when updating packages after using Pixi.

Fixed by disabling file locking with `conda config --set no_lock true`. 
[See here](https://github.com/conda/conda/issues/13534#issuecomment-2890328744).

