## Example MalariaGEN analysis: CNV https://malariagen.github.io/vector-data/ag3/examples/cnv-explore.html
## activate conda environment in the terminal below with installed python packages
## then use here for import

# %%

import malariagen_data
from bisect import bisect_left, bisect_right
import numpy as np
import dask.array as da
from dask.diagnostics import ProgressBar
import bokeh.io as bkio
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

# %% to make interactive plots with bokeh we need to set up plotting with boken 
bkio.output_notebook()

# %% setup access to Ag3 data in Google Cloud
ag3 = malariagen_data.Ag3()
ag3

# %%
