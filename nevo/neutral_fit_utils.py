import os
import numpy as np
from scipy.stats import beta 

def beta_cdf(p, N, m):
    '''Expected long term distribution under the 
        neutral model (truncated cumulative beta-distribution)
    
    Copyright
    ---------
    Github: https://github.com/misieber/neufit
    Theroy as described in [1]_. and [2]_.
    
    Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions 
    are met:

    1. Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above 
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    References
    ----------
    .. [1] Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. 
    and Curtis, T. P. (2006). Quantifying the roles of immigration 
    and chance in shaping prokaryote community structure. Environmental 
    Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x
    
    .. [2] Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., 
    Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, 
    S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. 
    (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243
    '''
    return beta.cdf(1.0, N*m*p, N*m*(1.0-p)) - beta.cdf(1.0/N, N*m*p, N*m*(1.0-p))

def subsample(counts, depth):
    '''Subsamples counts to uniform depth, dropping all samples 
        without enough depth 
    
    Copyright
    ---------
    Github: https://github.com/misieber/neufit
    Theroy as described in [1]_. and [2]_.
    
    Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions 
    are met:

    1. Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above 
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    References
    ----------
    .. [1] Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. 
    and Curtis, T. P. (2006). Quantifying the roles of immigration 
    and chance in shaping prokaryote community structure. Environmental 
    Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x
    
    .. [2] Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., 
    Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, 
    S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. 
    (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243    
    '''
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
    '''Argparser type: non-negative int
    
    Copyright
    ---------
    Github: https://github.com/misieber/neufit
    Theroy as described in [1]_. and [2]_.
    
    Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions 
    are met:

    1. Redistributions of source code must retain the above copyright 
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above 
       copyright notice, this list of conditions and the following 
       disclaimer in the documentation and/or other materials provided 
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
    "AS IS" ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    References
    ----------
    .. [1] Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. 
    and Curtis, T. P. (2006). Quantifying the roles of immigration 
    and chance in shaping prokaryote community structure. Environmental 
    Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x
    
    .. [2] Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., 
    Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, 
    S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. 
    (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243
    '''
    
    nnint = int(arg)
    if nnint < 0:
        raise ArgumentTypeError(arg + ' < 0, must be non-negative')
    return nnint