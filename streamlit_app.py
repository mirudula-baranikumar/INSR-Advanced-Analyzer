import streamlit as st
import requests, re, csv, io
from datetime import datetime

st.set_page_config(page_title="INSR Advanced Analyzer", page_icon="🧬", layout="wide")
ACC="P06213"
URL=f"https://rest.uniprot.org/uniprotkb/{ACC}.json"

st.title("🧬 INSR Advanced Analyzer")
st.caption("Motif and Domain Analysis of Insulin Receptor Variants Associated with Insulin Resistance")
st.info("Educational bioinformatics application. Results are not a clinical diagnosis or validated pathogenicity prediction.")

@st.cache_data(ttl=3600)
def get_uniprot():
    r=requests.get(URL,headers={"User-Agent":"INSR-Advanced-Analyzer-Web/1.0"},timeout=30)
    r.raise_for_status()
    d=r.json()
    seq=d.get("sequence",{}).get("value","")
    if not seq: raise ValueError("No sequence returned by UniProt.")
    return d,seq

def features_from(d):
    keep={"Domain","Region","Motif","Repeat","Chain","Topological domain","Transmembrane",
          "Signal peptide","Peptide","Active site","Binding site","Site","Zinc finger"}
    out=[]
    for f in d.get("features",[]):
        loc=f.get("location",{})
        start=loc.get("start",{}).get("value")
        end=loc.get("end",{}).get("value")
        if start is None or end is None or f.get("type") not in keep: continue
        out.append({"Type":f.get("type",""),"Annotation":f.get("description") or f.get("type",""),
                    "Start":int(start),"End":int(end)})
    return sorted(out,key=lambda x:x["Start"])

PATTERNS=[
    ("Protein kinase glycine-rich loop",r"[LIV]G.G..G","Simplified kinase-associated pattern."),
    ("Protein kinase HRD motif",r"HRD","Catalytic-loop signature."),
    ("Protein kinase DFG motif",r"DFG","Activation-segment signature."),
    ("Protein kinase APE motif",r"APE","Kinase activation-segment-associated signature.")
]
def scan(seq):
    rows=[]
    for name,pat,note in PATTERNS:
        for m in re.finditer(pat,seq):
            rows.append({"Motif":name,"Pattern":pat,"Start":m.start()+1,"End":m.end(),
                         "Match":m.group(),"Note":note})
    return rows

def parse_variant(s,seq):
    s=s.strip()
    m=re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])",s,re.I)
    if not m: raise ValueError("Enter a substitution like R1158W (one-letter amino-acid format).")
    ref,pos,alt=m.group(1).upper(),int(m.group(2)),m.group(3).upper()
    if pos<1 or pos>len(seq): raise ValueError(f"Position must be between 1 and {len(seq)}.")
    actual=seq[pos-1]
    if ref!=actual: raise ValueError(f"Reference mismatch at position {pos}: sequence contains {actual}, not {ref}. Check isoform/numbering.")
    if ref==alt: raise ValueError("Reference and alternate residues are identical.")
    return ref,pos,alt

with st.sidebar:
    st.header("Analysis settings")
    st.write("Reference protein: **Human insulin receptor**")
    st.write("UniProt accession: **P06213**")
    st.write("Data source: [UniProt REST API](https://rest.uniprot.org/uniprotkb/P06213.json)")
    if st.button("Refresh UniProt data"):
        get_uniprot.clear()
        st.rerun()

try:
    data,seq=get_uniprot()
    feats=features_from(data)
    motifs=scan(seq)
    protein_name=data.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value","Insulin receptor")
    organism=data.get("organism",{}).get("scientificName","Homo sapiens")
except Exception as e:
    st.error(f"Could not retrieve UniProt data: {e}")
    st.warning("Check the hosting environment's internet access and refresh the page.")
    st.stop()

c1,c2,c3=st.columns(3)
c1.metric("Sequence length",f"{len(seq):,} aa")
c2.metric("Annotated features",len(feats))
c3.metric("Scanned motif matches",len(motifs))
st.success(f"Loaded {protein_name} | {organism} | UniProt {ACC}")

tab1,tab2,tab3,tab4,tab5=st.tabs(["Overview","Motif analysis","Domain analysis","Variant analysis","Sequence & reports"])
with tab1:
    st.subheader("Protein overview")
    st.write("The canonical UniProt protein record is retrieved live. Sequence and feature annotations may change as UniProt updates its record.")
    st.json({"Protein":protein_name,"Organism":organism,"Accession":ACC,"Sequence length":len(seq),
             "Retrieved at":datetime.now().strftime("%Y-%m-%d %H:%M")})
    st.markdown("**Workflow**")
    st.markdown("1. Retrieve the reference insulin receptor sequence.\n2. Scan selected kinase-associated sequence signatures.\n3. Display UniProt-annotated protein features.\n4. Validate a substitution and map it to overlapping features/motifs.\n5. Export tables and reports.")
with tab2:
    st.subheader("Conserved motif-pattern scan")
    st.caption("These are simplified sequence-pattern searches, not a full PROSITE or ScanProsite implementation.")
    if motifs: st.dataframe(motifs,use_container_width=True,hide_index=True)
    else: st.warning("No matches were found for the selected patterns in this sequence.")
    import pandas as pd
    st.download_button("Download motif results CSV",pd.DataFrame(motifs).to_csv(index=False).encode("utf-8"),
                       "insr_motif_results.csv","text/csv",disabled=not motifs)
with tab3:
    st.subheader("UniProt annotated protein features")
    st.caption("The returned feature list may include domains, regions, sites and other annotations; not every feature is a domain.")
    if feats: st.dataframe(feats,use_container_width=True,hide_index=True)
    else: st.warning("No matching feature annotations were returned.")
    import pandas as pd
    st.download_button("Download feature annotations CSV",pd.DataFrame(feats).to_csv(index=False).encode("utf-8"),
                       "insr_feature_annotations.csv","text/csv",disabled=not feats)
    if feats:
        st.markdown("**Feature position overview**")
        st.caption("Start/end coordinates are based on the UniProt canonical sequence.")
        st.dataframe([{"Feature":x["Annotation"],"Region":f'{x["Start"]}–{x["End"]}',"Length":x["End"]-x["Start"]+1} for x in feats],
                     use_container_width=True,hide_index=True)
with tab4:
    st.subheader("Variant analysis")
    variant=st.text_input("Protein substitution",value="R1158W",placeholder="e.g. R1158W")
    if st.button("Analyze variant",type="primary"):
        try:
            ref,pos,alt=parse_variant(variant,seq)
            overlap=[x for x in feats if x["Start"]<=pos<=x["End"]]
            nearby=[x for x in feats if x["Start"]<=pos+10 and x["End"]>=pos-10 and x not in overlap]
            mo=[x for x in motifs if x["Start"]<=pos<=x["End"]]
            st.success(f"Validated substitution: {ref}{pos}{alt}")
            a,b,c=st.columns(3)
            a.metric("Reference residue",ref); b.metric("Position",pos); c.metric("Alternate residue",alt)
            st.write("**Features containing this position**")
            if overlap: st.dataframe(overlap,use_container_width=True,hide_index=True)
            else: st.write("No retrieved UniProt feature directly spans this position.")
            st.write("**Nearby features (within 10 residues)**")
            if nearby: st.dataframe(nearby,use_container_width=True,hide_index=True)
            else: st.write("No nearby annotated feature found.")
            st.write("**Scanned motifs containing this position**")
            if mo: st.dataframe(mo,use_container_width=True,hide_index=True)
            else: st.write("No scanned motif match directly spans this position.")
            st.warning("Feature overlap is not proof of functional disruption. This tool does not independently establish pathogenicity. Check ClinVar, UniProt evidence and primary literature.")
            st.session_state["last_variant"]={"Variant":f"{ref}{pos}{alt}","Position":pos,"Reference":ref,"Alternate":alt,
                "Overlapping features":overlap,"Nearby features":nearby,"Motif matches":mo}
        except Exception as e: st.error(str(e))
with tab5:
    st.subheader("Reference FASTA")
    fasta=f">{ACC}|canonical|{organism}|length={len(seq)}\n"+"\n".join(seq[i:i+60] for i in range(0,len(seq),60))
    st.download_button("Download reference FASTA",fasta,"INSR_P06213.fasta","text/plain")
    st.text_area("FASTA sequence",fasta,height=220)
    report=["INSR ADVANCED ANALYZER",f"UniProt accession: {ACC}",f"Protein: {protein_name}",
            f"Organism: {organism}",f"Sequence length: {len(seq)}",f"Generated: {datetime.now().isoformat(timespec='minutes')}","","MOTIF MATCHES"]
    report += [f'{x["Motif"]} | {x["Pattern"]} | {x["Start"]}-{x["End"]} | {x["Match"]}' for x in motifs]
    report += ["","UNIPROT FEATURES"]
    report += [f'{x["Type"]} | {x["Annotation"]} | {x["Start"]}-{x["End"]}' for x in feats]
    if "last_variant" in st.session_state:
        report += ["","LATEST VARIANT",str(st.session_state["last_variant"])]
    report += ["","Limitations: educational research use only; motif patterns are simplified; feature overlap is not a pathogenicity prediction."]
    st.download_button("Download complete TXT report","\n".join(report),"INSR_analysis_report.txt","text/plain")

st.divider()
st.caption("Sources: UniProtKB P06213. Motif pattern scans are simplified educational searches. Verify all findings against authoritative annotations and publications.")
