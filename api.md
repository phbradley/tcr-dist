# API branch
---


# Single TCR example
Here's an example of how you can process a single nucleotide sequence:

```python
betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'

betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'

import tcrdist as td
chain = td.processing.processNT('human', 'B', betaNT, betaQuals)
```

# Pipeline example
Here's an example of a TCR pipeline with `tcrdist`:

```python
import tcrdist as td
psDf = td.datasets.loadPSData('human_pairseqs_v1')
probDf = td.processing.computeProbs(psDf)
psDf = psDf.join(probDf)
clonesDf = td.processing.identifyClones(psDf)
```

# Testing
To run all tests:

`py.test tcrdist`

or to pass all print statements through for debugging:

`py.test tcrdist -s`

Note that importing and running functions in `tcrdist` generates a log file, `tcrdist.log`. You can set the logging `level` by modifying the parameter in the base `__init__.py` file using one of the following:

```
logging.ERROR
logging.WARNING
logging.INFO
logging.DEBUG
```

For details see documentation on the `logging` [module](https://docs.python.org/2/library/logging.html).
