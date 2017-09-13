
def parse_tsv_line(line,infields):
    if line[-1]=='\n': line = line[:-1] #doh
    l = line.split('\t')
    assert len(l) == len(infields)
    vals = {}
    for tag,val in zip(infields,l):
        vals[tag] = val
    return vals

def make_tsv_line(vals,outfields,empty_string_replacement=''):
    """Does not have the \n at the end"""
    l = []
    for tag in outfields:
        val = vals[tag]
        if type(val) is str:
            if empty_string_replacement and not val:
                l.append( empty_string_replacement )
            else:
                l.append(val)
        else:
            l.append(str(val))
    return '\t'.join( l )



def parse_tsv_file( filename, key_fields=[], store_fields=[], save_l=False ):
    if not key_fields and not store_fields:
        save_l = True
    D = {}
    L = []

    infields = []
    for line in open( filename,'rU'):
        if not infields:
            if line[0] == '#':
                infields = line[1:-1].split('\t')
            else:
                infields = line[:-1].split('\t')
            continue
        assert infields

        l = parse_tsv_line( line[:-1], infields )

        if store_fields:
            dats = [ l[x] for x in store_fields ]
            if save_l:
                dats.append( l )
        else:
            assert save_l
            dats = l


        if key_fields:
            subd = D
            for k in key_fields[:-1]:
                tag = l[k]
                if tag not in subd: subd[tag] = {}
                subd = subd[tag]

            final_tag = l[ key_fields[-1] ]
            if final_tag not in subd: subd[final_tag] = []

            subd[final_tag].append( dats )
        else:
            L.append( dats )

    if key_fields:
        return D
    else:
        return L

