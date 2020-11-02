def Collapse(x,y):

    import numpy as np
    from pdb import set_trace

    x_new = np.unique(x)
    y_new = np.array([])
    for value in x_new:
        indices = x == value
        y_new_value = sum(y[indices])
        y_new = np.append(y_new,y_new_value)

    return(x_new,y_new)
