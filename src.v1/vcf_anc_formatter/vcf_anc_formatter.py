import sys, gzip, re, os
import pprint as pp

if len(sys.argv) != 2 :
    exit("usage: zcat [input: Flare's output vcf.gz] | python3 vcf_anc_formatter.py [output prefix]")


outpref = sys.argv[1]


regex = re.compile(r'ANCESTRY=<(.*?)>')

ancs = {}
output = {}

headline = ''

for line in sys.stdin :
    if not line: break
    ## read ancestry code
    if line[0:2] == "##" :
        if "fileformat" in line : headline += line
        if "source" in line: headline += line
        if "ID=GT" in line: headline += line

        ret = regex.search(line)
        if ret :
            anc_str = ret.group(1)
            for x in anc_str.split(',') :
                anc,code = x.split('=')
                ancs[code] = anc
            for code in ancs :
                output[code] = {}
                output[code]['genofile'] = outpref + ".anc" + code + ".vcf"
                output[code]['hapfile'] = outpref + ".hapanc" + code + ".vcf"
                output[code]['ancestry'] = ancs[code]
            pp.pprint(output)
    elif "#CHROM" in line :
        headline += line
        if not ancs :
            exit("wrong input, no ancestry code information")
        for code in ancs:
            with open(output[code]['genofile'], 'w') as fp: fp.write(headline)
            with open(output[code]['hapfile'], 'w') as fp: fp.write(headline)
    else :
        # read in data
        llst = line.split()
        meta_col = llst[0:8]
        meta_col.append('GT')
        genos = {}
        haps = {}
        for code in ancs :
            genos[code] = []
            haps[code] = []
        for x in llst[9:] :
            y = x.split(':')
            g = y[0]
            anc_codes = y[1:]
            for code in ancs:
                gx = ''
                hx = ''
                if anc_codes[0] == code :
                    gx += g[0]
                    hx += '1'
                else :
                    gx += '.'
                    hx += '0'
                gx += g[1]
                hx += g[1]
                if anc_codes[1] == code :
                    gx += g[2]
                    hx += '1'
                else :
                    gx += '.'
                    hx += '0'
                genos[code].append(gx)
                haps[code].append(hx)
        #output
        for code in ancs:
            ogeno = '\t'.join(meta_col) + '\t' + '\t'.join(genos[code]) + '\n'
            ohap = '\t'.join(meta_col) + '\t' + '\t'.join(haps[code]) + '\n'
            with open(output[code]['genofile'], 'at') as fp: fp.write(ogeno)
            with open(output[code]['hapfile'], 'at') as fp: fp.write(ohap)


# compressed the file
for code in ancs:
    with open(output[code]['genofile'], 'rb') as f_in, gzip.open(output[code]['genofile']+'.gz', 'wb') as f_out: f_out.writelines(f_in)
    with open(output[code]['hapfile'], 'rb') as f_in, gzip.open(output[code]['hapfile']+'.gz', 'wb') as f_out: f_out.writelines(f_in)
    os.remove(output[code]['genofile']) 
    os.remove(output[code]['hapfile']) 
