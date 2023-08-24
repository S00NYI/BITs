import itertools

def barcodes_N(barcode, n):
    positions = range(len(barcode))
    for pos in itertools.combinations(positions, n):
        new_barcode = list(barcode)
        for p in pos:
            new_barcode[p] = 'N'
        yield ''.join(new_barcode)


AM025_Input = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC']
JL0361_Input = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT']
JL0361_Enrich = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT']
JL0380_Fraction = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA']
JL0388_Input = ['ATCACG', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTCCGC', 'GTGAAA']
JL0388_Enrich = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA', 'GTCCGC', 'GTGAAA']
JL1024_Pool = ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA', 'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTCCGC', 'GTGAAA']

Ns_allowed = 2

barcodes = JL1024_Pool
new_barcodes = []
for barcode in barcodes:
    for combination in barcodes_N(barcode, Ns_allowed):
        new_barcodes.append(combination)

## Check whether all the barcodes with Ns are unique or not.
## If True, there is no overlap.
## IF False, there is overlap (we lose uniqueness).
len(new_barcodes) == len(list(set(new_barcodes)))
