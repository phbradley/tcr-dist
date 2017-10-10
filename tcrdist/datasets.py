import os.path as op
import logging
import inspect

from . import processing

logger = logging.getLogger('datasets.py')

def loadPSData(name=None):
    """Used to load paired sequence data without needing to be aware of the full path.
    Call without parameters to see a list of datasets available."""
    
    datasets = dict(human_pairseqs_v1='human_pairseqs_v1.tsv',
                    mouse_pairseqs_v1='mouse_pairseqs_v1.tsv',
                    test_human_pairseqs='test_human_pairseqs.tsv',
                    test_mouse_pairseqs='test_mouse_pairseqs.tsv')
    
    if name is None or name == '':
        print('Available datsets:')
        for k in datasets.keys():
            print('  %s' % k)
        return
    
    try:
        fn = datasets[name]
    except KeyError:
        logger.error('Could not find dataset: %s' % name)
    
    if fn[-3:] == 'tsv':
        delimiter = '\t'
    else:
        delimiter = ','
    
    datasetsPath = op.join(op.split(inspect.stack()[0][1])[0], 'datasets')
    filename = op.join(datasetsPath, fn)
    
    print('Reading from %s' % filename)
    
    if 'mouse' in fn:
        organism = 'mouse'
    else:
        organism = 'human'
    
    psDf = processing.readPairedSequences(organism, filename)
    return psDf
    
    
        