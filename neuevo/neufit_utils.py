# neufit: Fit a neutral community model to species abundances, e.g. from an OTU table
#
# For the theory behind this see Sloan et al, Environ Microbiol 2006 8:732-740.
# To run on the example simulation data: python neufit.py sim_data.csv
# To link with the mock taxonomy use the -t sim_taxonomy.csv option
#
# Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Github: https://github.com/misieber/neufit

import os
import numpy as np
from scipy.stats import beta 
from datetime import datetime

def beta_cdf(p, N, m):
    # Expected long term distribution under the neutral model (truncated cumulative beta-distribution)
    return beta.cdf(1.0, N*m*p, N*m*(1.0-p)) - beta.cdf(1.0/N, N*m*p, N*m*(1.0-p))

def subsample(counts, depth):
    # Subsamples counts to uniform depth, dropping all samples without enough depth
    for sample in counts:
        if counts[sample].sum() >= depth:
            flattened = np.repeat(np.arange(counts[sample].size), counts[sample])
            subsample = np.random.choice(flattened, depth, replace=False)
            counts[sample] = np.bincount(subsample, minlength=counts[sample].size)
        else:
            #CG: changed the following print statment from Python2 to Python3
            print('dropping sample ' + sample + ' with ' + str(counts[sample].sum()) + ' reads < ' + str(depth))
            counts = counts.drop(sample, axis=1)
    return counts

def non_negative_int(arg):
    # Argparser type: non-negative int
    nnint = int(arg)
    if nnint < 0:
        raise ArgumentTypeError(arg + ' < 0, must be non-negative')
    return nnint