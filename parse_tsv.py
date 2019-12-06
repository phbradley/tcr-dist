
def parse_tsv_line(line,infields,sep='\t'):
    if line[-1]=='\n': line = line[:-1] #doh
    l = line.split(sep)
    assert len(l) == len(infields)
    vals = {}
    for tag,val in zip(infields,l):
        vals[tag] = val
    return vals

def make_tsv_line(vals,outfields,empty_string_replacement='',sep='\t'):
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
    return sep.join( l )



def parse_tsv_file( filename, key_fields=[], store_fields=[], save_l=False, sep='\t' ):
    if not key_fields and not store_fields:
        save_l = True
    D = {}
    L = []

    infields = []
    for line in open( filename,'rU'):
        if not infields:
            if line[0] == '#':
                infields = line[1:-1].split(sep)
            else:
                infields = line[:-1].split(sep)
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

# silly
def safely_split_csv_line( line, text_encapsulator='"' ):
    newcomma = "COMMA!!DUDE!!" # any string that's not present in line will work
    assert newcomma not in line
    l = list(line)
    assert len(l) == len(line)
    in_quote = False
    newl = l[:]
    for i,a in enumerate(l):
        if a == text_encapsulator:
            in_quote = not in_quote
        else:
            if in_quote and a == ',':
                newl[i] = newcomma
    assert not in_quote
    #print ''.join(newl)

    l = ( ''.join( newl ) ).split(',')
    newl = l[:]
    for i,a in enumerate(l):
        if newcomma in a:
            newl[i] = l[i].replace(newcomma,',')

    return newl
