#!/usr/bin/env python3
import sys, gzip, argparse
from collections import defaultdict

def openz(p): return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')

def parse_attrs(s):
    d={}
    for it in s.strip().split(';'):
        it=it.strip()
        if not it: continue
        if ' ' in it:
            k,v=it.split(' ',1); v=v.strip().strip('"')
        elif '=' in it:
            k,v=it.split('=',1); v=v.strip().strip('"')
        else:
            continue
        d[k]=v
    return d

def merge(iv):
    if not iv: return []
    iv=sorted(iv)
    out=[]; s,e=iv[0]
    for a,b in iv[1:]:
        if a<=e:
            if b>e: e=b
        else:
            out.append((s,e)); s,e=a,b
    out.append((s,e))
    return out

def utr5_from_exon_cds(exons, cds, strand):
    if not exons or not cds: return []
    exons=merge(exons); cds=merge(cds)
    if strand=='+':
        tss=min(s for s,_ in exons)
        cds_start=min(s for s,_ in cds)
        if cds_start<=tss: return []
        window=[(tss, cds_start)]
    else:
        tss=max(e for _,e in exons)
        cds_end=max(e for _,e in cds)
        if cds_end>=tss: return []
        window=[(cds_end, tss)]
    out=[]; i=j=0
    while i<len(exons) and j<len(window):
        a1,a2=exons[i]; b1,b2=window[j]
        s=max(a1,b1); e=min(a2,b2)
        if s<e: out.append((s,e))
        if a2<=b2: i+=1
        else: j+=1
    return merge(out)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gene-name-attr", dest="gene_name_attr", default="gene_name",
                    help="GTF attribute to use for gene label; falls back to gene_id if missing")
    ap.add_argument("--out-bed", required=True, help="BED6 of union 5'UTR per gene")
    args=ap.parse_args()

    tx_exon=defaultdict(list)
    tx_cds =defaultdict(list)
    tx_gene={}
    with openz(args.gtf) as fh:
        for line in fh:
            if not line or line[0]=='#': continue
            f=line.rstrip('\n').split('\t')
            if len(f)<9: continue
            chrom,_,feat,start,end,_,strand,_,attr=f
            if feat not in ("exon","CDS","transcript"): continue
            a=parse_attrs(attr)
            tx=a.get("transcript_id") or a.get("transcript") or a.get("transcriptId")
            if not tx: continue
            # prefer gene_name_attr, then gene_name, then gene_id
            gene = a.get(args.gene_name_attr) or a.get("gene_name") or a.get("gene_id") or a.get("gene")
            if not gene: continue
            if feat=="transcript":
                tx_gene[tx]=gene
                continue
            s=int(start)-1; e=int(end)
            key=(chrom,strand,tx)
            if feat=="exon":
                tx_exon[key].append((s,e)); tx_gene[tx]=gene
            elif feat=="CDS":
                tx_cds[key].append((s,e)); tx_gene[tx]=gene

    from collections import defaultdict as dd
    per_gene_cs=dd(lambda: dd(list))
    for (chrom,strand,tx), ex in tx_exon.items():
        cds=tx_cds.get((chrom,strand,tx))
        if not cds: continue
        utr=utr5_from_exon_cds(ex, cds, strand)
        if utr:
            g=tx_gene.get(tx)
            if g: per_gene_cs[g][(chrom,strand)].extend(utr)

    out=[]
    for g, bycs in per_gene_cs.items():
        best=None; best_len=-1; best_cs=None
        for cs, iv in bycs.items():
            m=merge(iv); L=sum(e-s for s,e in m)
            if L>best_len: best_len=L; best=m; best_cs=cs
        if best and best_len>0:
            chrom,strand=best_cs
            for s,e in best:
                out.append((chrom,s,e,g,0,strand))

    out.sort(key=lambda x:(x[0],x[1],x[2],x[3]))
    with open(args.out_bed,'w') as oh:
        for chrom,s,e,g,score,strand in out:
            oh.write(f"{chrom}\t{s}\t{e}\t{g}\t0\t{strand}\n")

if __name__=="__main__":
    main()

