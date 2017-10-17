import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('objects.py')

from . import translation
from .blosum import bsd4,blosum

class DotDict(dict):
    def __getattr__(self, item):
        if item in self:
            return self[item]
        logger.error('No attribute %s!' % item)
        raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value
    def __str__(self):
        return pd.Series(self).to_string()
    def to_series(self):
        return pd.Series(self)

class TCRClone(DotDict):
    """Object that contains all info for a single TCR clone.
    As a DotDict attributes can be accessed using dot notation
    or standard dict key access."""

    alphaAA = ''
    betaAA = ''
    gammaAA = ''
    deltaAA = ''

    
    subjid = ''
    cloneid = ''
    epitope = ''

    def __init__(self, chain1, chain2, **kwargs):
        for k in chain1.keys():
            self[k] = chain1[k]
        for k in chain2.keys():
            self[k] = chain2[k]

        for k in kwargs.keys():
            self[k] = kwargs[k]

class TCRChain(DotDict):
    def __init__(self, **kwargs):
        for k in kwargs.keys():
            self[k] = kwargs[k]

class TCR_Gene:
    cdrs_sep = ';'
    gap_character = '.'

    def __init__( self, l):
        self.id = l['id']
        self.organism = l['organism']
        self.chain = l['chain']
        self.region = l['region']
        self.nucseq = l['nucseq']
        self.alseq = l['aligned_protseq']
        if pd.isnull(l['cdrs']):
            self.cdrs = []
            self.cdr_columns = []
        else:
            self.cdrs = l['cdrs'].split(self.cdrs_sep)
            ## these are still 1-indexed !!!!!!!!!!!!!!
            self.cdr_columns = [ map( int, x.split('-')) for x in l['cdr_columns'].split(self.cdrs_sep) ]
        frame = l['frame']
        assert frame in [1, 2, 3]
        self.nucseq_offset = frame - 1 ## 0, 1 or 2 (0-indexed for python)
        self.protseq = translation.get_translation( self.nucseq, frame )[0]
        assert self.protseq == self.alseq.replace(self.gap_character,'')
        # sanity check
        if self.cdrs:
            assert self.cdrs == [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]

class DistanceParams(DotDict):
    def __init__(self, config_string=None):
        self.gap_penalty_v_region = 4
        self.gap_penalty_cdr3_region = 12 # same as gap_penalty_v_region=4 since weight_cdr3_region=3 is not applied
        self.weight_v_region = 1
        self.weight_cdr3_region = 3
        self.distance_matrix = bsd4
        self.align_cdr3s = False
        self.trim_cdr3s = True
        self.scale_factor = 1.0
        if config_string:
            l = config_string.split(',')
            for tag,val in [x.split(':') for x in l ]:
                if tag == 'gap_penalty_cdr3_region':
                    self.gap_penalty_cdr3_region = float(val)
                elif tag == 'gap_penalty_v_region':
                    self.gap_penalty_v_region = float(val)
                elif tag == 'weight_cdr3_region':
                    self.weight_cdr3_region = float(val)
                elif tag == 'weight_v_region':
                    self.weight_v_region = float(val)
                elif tag == 'scale_factor':
                    self.scale_factor = float(val)
                elif tag == 'align_cdr3s':
                    assert val in ['True','False']
                    self.align_cdr3s = ( val == 'True' )
                elif tag == 'trim_cdr3s':
                    assert val in ['True','False']
                    self.trim_cdr3s = ( val == 'True' )
                else:
                    logger.error('unrecognized tag: %s' % tag)
                    assert False
            logger.info('config_string: {}\nself:\n{}'.format(config_string, self))
