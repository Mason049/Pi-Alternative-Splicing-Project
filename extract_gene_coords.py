#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, gzip, os, sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description=(
        "Extract genomic coordinates by matching each row's 'seqname' in the associated GFF file. "
        "GFFs are located by exact basename under --gff-root (recursive). Optionally require co-location with genome."
    ))
    p.add_argument("--metadata", required=True)
    p.add_argument("--gff-root", required=True, help="Root directory containing GFF files (searched recursively).")
    p.add_argument("--pairing", choices=["basename","co-located"], default="basename",
                   help="basename: accept any GFF path by basename; co-located: require same dir also has the genome basename (or twin).")
    p.add_argument("--gff-col", default="GFF File Name")
    p.add_argument("--genome-col", default="Genome File Name")
    p.add_argument("--cds-col", default="CDS File Name")
    p.add_argument("--proteome-col", default="Proteome File Name")
    p.add_argument("--seqname-col", default="seqname")
    p.add_argument("--id-col", default="ID")
    p.add_argument("--delimiter", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--feature-priority", default="mRNA,gene,transcript,CDS,exon")
    p.add_argument("--match-mode", choices=["strict-id","id+name"], default="strict-id")
    return p.parse_args()

def open_possible_gzip(path:str):
    return gzip.open(path,'rt') if path.endswith('.gz') else open(path,'r')

def parse_gff_attributes(attr_str:str)->Dict[str,str]:
    attrs={}
    if not attr_str or attr_str=='.': return attrs
    for field in attr_str.strip().split(';'):
        if not field: continue
        if '=' in field: key,val=field.split('=',1)
        elif ' ' in field:
            parts=field.split(' ',1); key,val=(parts[0],parts[1]) if len(parts)==2 else (field,'')
        else: key,val=field,''
        attrs[key.strip()]=val.strip()
    return attrs

class GFFIndex:
    def __init__(self,gff_path:str):
        self.gff_path=gff_path
        self.id_to_rec:Dict[str,Tuple[str,str,str,int,int,str,Dict[str,str]]]={}
        self.name_to_ids:Dict[str,List[str]]=defaultdict(list)
        self.child_to_parents:Dict[str,List[str]]=defaultdict(list)
        self.type_by_id:Dict[str,str]={}
    def build(self):
        with open_possible_gzip(self.gff_path) as fh:
            for ln in fh:
                if not ln or ln.startswith('#'): continue
                parts=ln.rstrip('\n').split('\t')
                if len(parts)<9: continue
                seqid,source,ftype,start,end,score,strand,phase=parts[:8]
                attrs_raw='\t'.join(parts[8:])
                try: start_i=int(start); end_i=int(end)
                except: continue
                attrs=parse_gff_attributes(attrs_raw)
                rec_id=attrs.get('ID')
                if rec_id:
                    self.id_to_rec[rec_id]=(seqid,source,ftype,start_i,end_i,strand,attrs)
                    self.type_by_id[rec_id]=ftype
                name=attrs.get('Name') or attrs.get('gene') or attrs.get('gene_name') or attrs.get('locus_tag')
                if name and rec_id: self.name_to_ids[name].append(rec_id)
                parents=attrs.get('Parent')
                if parents and rec_id:
                    for p in parents.split(','): self.child_to_parents[rec_id].append(p)
    def get_record(self,rec_id:str): return self.id_to_rec.get(rec_id)
    def get_parents(self,rec_id:str)->List[str]: return self.child_to_parents.get(rec_id,[])
    def get_type(self,rec_id:str)->Optional[str]: return self.type_by_id.get(rec_id)
    def lookup_by_name(self,name:str)->List[str]: return self.name_to_ids.get(name,[])
    def ascend_to_gene(self,rec_id:str)->Tuple[str,Optional[Tuple[str,str,str,int,int,str,Dict[str,str]]]]:
        visited=set(); cur=rec_id
        while cur and cur not in visited:
            visited.add(cur); rec=self.get_record(cur)
            if not rec: break
            _,_,ftype,*_=rec
            if ftype=='gene': return cur,rec
            parents=self.get_parents(cur)
            if not parents: return cur,rec
            gene_parents=[p for p in parents if self.get_type(p)=='gene']
            cur=gene_parents[0] if gene_parents else parents[0]
        return rec_id,self.get_record(rec_id)

def choose_best_id(cands:List[str], index:GFFIndex, prio_list:List[str])->Optional[str]:
    if not cands: return None
    prio={t:i for i,t in enumerate(prio_list)}; scored=[]
    for cid in cands:
        ftype=index.get_type(cid) or ''
        scored.append((prio.get(ftype,len(prio)+1), len(cid), cid))
    scored.sort(); return scored[0][2]

def find_match_for_label(label:str, index:GFFIndex, prio_list:List[str], mode:str)->Tuple[Optional[str],str]:
    if index.get_record(label): return label,'ID_exact'
    base=label
    if '.' in label:
        base=label.rsplit('.',1)[0]
        if index.get_record(base): return base,'ID_base'
    if mode=='id+name':
        ids=index.lookup_by_name(label)
        if ids: return choose_best_id(ids,index,prio_list),'Name_exact'
        if base!=label:
            ids2=index.lookup_by_name(base)
            if ids2: return choose_best_id(ids2,index,prio_list),'Name_base'
    return None,'no_match'

def twin_name(basename:str)->List[str]:
    return [basename[:-3]] if basename.endswith('.gz') else [basename+'.gz']

def find_gff_by_basename(root:str, gff_name:str, genome_name:Optional[str], pairing:str)->Optional[str]:
    cands=[]
    for dirpath,_,files in os.walk(root):
        if gff_name in files:
            if pairing=='co-located' and genome_name:
                twins=twin_name(genome_name)
                if (genome_name in files) or any(t in files for t in twins):
                    cands.append(os.path.join(dirpath,gff_name))
            else:
                cands.append(os.path.join(dirpath,gff_name))
    if not cands: return None
    cands.sort(key=lambda p:(p.count(os.sep), len(p), p))
    return cands[0]

def main():
    args=parse_args()
    read_kwargs={}
    if args.delimiter: read_kwargs['sep']=args.delimiter
    try: meta=pd.read_csv(args.metadata,**read_kwargs)
    except Exception: meta=pd.read_csv(args.metadata, sep='\t')
    for c in [args.seqname_col, args.id_col, args.gff_col]:
        if c not in meta.columns:
            sys.stderr.write(f"ERROR: required column '{c}' not in metadata\n"); sys.exit(2)

    meta['_gff_name']=meta[args.gff_col].astype(str).apply(os.path.basename)
    meta['_genome_name']=meta[args.genome_col].astype(str).apply(os.path.basename) if args.genome_col in meta.columns else ''
    meta['_gff_path']=''
    seen_warn=set()
    for i,r in meta.iterrows():
        gff_name=r['_gff_name']
        genome_name=r['_genome_name'] if pd.notna(r.get('_genome_name','')) and r.get('_genome_name','')!='' else None
        path=find_gff_by_basename(args.gff_root, gff_name, genome_name, args.pairing)
        if not path:
            key=(gff_name, genome_name, args.pairing)
            if key not in seen_warn:
                sys.stderr.write(f"WARN: Could not find GFF '{gff_name}' (pairing={args.pairing}; genome='{genome_name}') under {args.gff_root}\n")
                seen_warn.add(key)
            meta.at[i,'_gff_path']=''
        else:
            meta.at[i,'_gff_path']=path

    groups=meta.groupby('_gff_path', dropna=False)
    prio=[t.strip() for t in args.feature_priority.split(',') if t.strip()]
    out_rows=[]
    for gff_path, group in groups:
        if not gff_path or not os.path.exists(gff_path):
            for _,r in group.iterrows():
                genome_val=str(r.get(args.genome_col,'')) if args.genome_col in r and pd.notna(r[args.genome_col]) else ''
                cds_val=str(r.get(args.cds_col,'')) if args.cds_col in r and pd.notna(r[args.cds_col]) else ''
                prot_val=str(r.get(args.proteome_col,'')) if args.proteome_col in r and pd.notna(r[args.proteome_col]) else ''
                out_rows.append({args.id_col:r[args.id_col], args.seqname_col:r[args.seqname_col], args.gff_col:r[args.gff_col],
                                 args.genome_col:genome_val, args.cds_col:cds_val, args.proteome_col:prot_val,
                                 'gff_path':'', 'match_how':'gff_missing', 'matched_feature_type':'','matched_feature_id':'',
                                 'gene_id':'','chrom':'','start':'','end':'','strand':''})
            continue

        index=GFFIndex(gff_path); index.build()
        for _,row in group.iterrows():
            label=str(row[args.seqname_col])
            genome_val=str(row.get(args.genome_col,'')) if args.genome_col in row and pd.notna(row[args.genome_col]) else ''
            cds_val=str(row.get(args.cds_col,'')) if args.cds_col in row and pd.notna(row[args.cds_col]) else ''
            prot_val=str(row.get(args.proteome_col,'')) if args.proteome_col in row and pd.notna(row[args.proteome_col]) else ''
            mid,how=find_match_for_label(label,index,prio,args.match_mode)
            if not mid:
                out_rows.append({args.id_col:row[args.id_col], args.seqname_col:label, args.gff_col:row[args.gff_col],
                                 args.genome_col:genome_val, args.cds_col:cds_val, args.proteome_col:prot_val,
                                 'gff_path':gff_path, 'match_how':how, 'matched_feature_type':'','matched_feature_id':'',
                                 'gene_id':'','chrom':'','start':'','end':'','strand':''}); continue
            gene_id, gene_rec=index.ascend_to_gene(mid)
            if gene_rec:
                chrom,source,ftype,start,end,strand,attrs=gene_rec
                mrec=index.get_record(mid); mtype=mrec[2] if mrec else ''
                out_rows.append({args.id_col:row[args.id_col], args.seqname_col:label, args.gff_col:row[args.gff_col],
                                 args.genome_col:genome_val, args.cds_col:cds_val, args.proteome_col:prot_val,
                                 'gff_path':gff_path, 'match_how':how, 'matched_feature_type':mtype, 'matched_feature_id':mid,
                                 'gene_id':gene_id, 'chrom':chrom, 'start':start, 'end':end, 'strand':strand})
            else:
                mrec=index.get_record(mid)
                if mrec:
                    chrom,source,ftype,start,end,strand,attrs=mrec
                    out_rows.append({args.id_col:row[args.id_col], args.seqname_col:label, args.gff_col:row[args.gff_col],
                                     args.genome_col:genome_val, args.cds_col:cds_val, args.proteome_col:prot_val,
                                     'gff_path':gff_path, 'match_how':how, 'matched_feature_type':ftype, 'matched_feature_id':mid,
                                     'gene_id':'', 'chrom':chrom, 'start':start, 'end':end, 'strand':strand})
                else:
                    out_rows.append({args.id_col:row[args.id_col], args.seqname_col:label, args.gff_col:row[args.gff_col],
                                     args.genome_col:genome_val, args.cds_col:cds_val, args.proteome_col:prot_val,
                                     'gff_path':gff_path, 'match_how':'unexpected_no_record', 'matched_feature_type':'','matched_feature_id':'',
                                     'gene_id':'','chrom':'','start':'','end':'','strand':''})

    out_df=pd.DataFrame(out_rows)
    cols=[args.id_col, args.seqname_col, args.gff_col, args.genome_col, args.cds_col, args.proteome_col,
          'gff_path','match_how','matched_feature_type','matched_feature_id','gene_id','chrom','start','end','strand']
    for c in cols:
        if c not in out_df.columns: out_df[c]=''
    out_df=out_df[cols]
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows to {args.out}", file=sys.stderr)

if __name__=='__main__':
    main()
