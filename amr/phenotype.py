"""Phenotype orchestration."""
from __future__ import annotations

from .engine import (
    AD_CAT,
    AD_MCF,
    AD_MIC,
    ANTIBIOTIC_ORDER,
    CATEGORIES,
    Counter,
    DATE_FIELDS,
    EST_CAT,
    EST_MCF,
    MIC_FIELDS,
    QC,
    ROUNDED,
    WGS_CAT,
    agreement_tables,
    build_mic_illustration_pptx,
    category,
    defaultdict,
    math,
    measurement,
    method_profiles,
    np,
    openpyxl,
    plt,
    style_workbook,
)


def phenotype_outputs(data, assemblies, root, audit):
    wb=openpyxl.Workbook()
    wide=wb.active
    wide.title="Isolate MICs"
    wide.append(["Canonical ID","Workbook IDs","Cohort"]+[f"{d} | {f} (µg/mL)" for d in ANTIBIOTIC_ORDER for f in MIC_FIELDS])
    details=wb.create_sheet("Method details")
    fields=(*MIC_FIELDS,*CATEGORIES,EST_MCF,AD_MCF,*DATE_FIELDS)
    details.append(["Canonical ID","Antibiotic",*[f+" (µg/mL)" if f in MIC_FIELDS else f for f in fields]])
    ids=wb.create_sheet("ID reconciliation")
    ids.append(["Canonical ID","Workbook IDs","Annotated source files","Assembly source files","Availability"])
    qc=wb.create_sheet("QC")
    qc.append(["Antibiotic",*[f+" (µg/mL)" for f in MIC_FIELDS],EST_CAT,AD_CAT,EST_MCF,AD_MCF])
    for key in sorted(set(assemblies)|set(data.ids)):
        raw=" | ".join(sorted(data.ids.get(key,[])))
        entry=assemblies.get(key,{})
        ids.append([key,raw," | ".join(p.name for p in entry.get("annotated_files",[])),
                    " | ".join(p.name for p in entry.get("plane_files",[])),
                    "both" if raw and entry else "phenotype only" if raw else "sequence only"])
        if not raw: continue
        wide.append([key,raw,"QC" if key==QC else "clinical"]+[data.values.get((key,d),{}).get(f) for d in ANTIBIOTIC_ORDER for f in MIC_FIELDS])
        for drug in ANTIBIOTIC_ORDER:
            row=data.values.get((key,drug),{})
            details.append([key,drug,*[row.get(f) for f in fields]])
            if key==QC: qc.append([drug,*[row.get(f) for f in (*MIC_FIELDS,EST_CAT,AD_CAT,EST_MCF,AD_MCF)]])
    trace=wb.create_sheet("Source cells")
    trace.append(["Canonical ID","Workbook ID","Antibiotic","Field","Source sheet","Cell","Source value","Source unit","Status"])
    for row in data.trace:trace.append(row)
    figdir=root/"MIC_Figures"
    figdir.mkdir(exist_ok=True)
    stats=wb.create_sheet("Paired agreement")
    stats.append(["Antibiotic","Exact quantitative pairs","Within one doubling dilution","Essential agreement (%)","Category pairs","Identical categories","Categorical agreement (%)","Log2 bias EST minus AD","Wilcoxon paired p"])
    concord=wb.create_sheet("Recorded WGS agreement")
    concord.append(["Antibiotic","Phenotype method","Paired clinical categories","Concordant","Agreement (%)"])
    for method in (ROUNDED, AD_MIC):
        counts=defaultdict(Counter)
        cat=EST_CAT if method==ROUNDED else AD_CAT
        for (key,d),row in data.values.items():
            if key==QC or row.get(method) is None: continue
            counts[d][str(row[method])]+=1
        if not counts:
            audit.add(method,"skipped","no clinical MICs in the source field")
            continue
        fig,axes=plt.subplots(4,2,figsize=(13,13),layout="constrained")
        for d,ax in zip(ANTIBIOTIC_ORDER,axes.flat):
            counts_d=counts[d]
            labels=sorted(counts_d,key=lambda v:(measurement(v)[1],v))
            ax.bar(range(len(labels)),[counts_d[k] for k in labels],color="#237c79")
            ax.set_xticks(range(len(labels)),labels,rotation=40,ha="right")
            ax.set_title(f"{d} (n={sum(counts_d.values())})")
            ax.set_ylabel("Clinical isolates")
            pairs=[(category(row.get(cat)),category(row.get(WGS_CAT))) for (k,drug),row in data.values.items()
                   if k!=QC and drug==d and category(row.get(cat)) and category(row.get(WGS_CAT))]
            agree=sum(a==b for a,b in pairs)
            concord.append([d,method,len(pairs),agree,100*agree/len(pairs) if pairs else None])
        fig.suptitle(method+" (µg/mL) • reported bounds retained")
        fig.savefig(figdir/("EST_rounded_MIC.png" if method==ROUNDED else "AD_MIC.png"),dpi=160)
        plt.close(fig)
    plotted=[]
    for drug in ANTIBIOTIC_ORDER:
        vals=[row for (key,d),row in data.values.items() if key!=QC and d==drug]
        pairs=[]
        cats=[]
        for row in vals:
            op1,a=measurement(row.get(ROUNDED));op2,b=measurement(row.get(AD_MIC))
            if a and b and op1 in ("","=") and op2 in ("","="):pairs.append((a,b))
            c1,c2=category(row.get(EST_CAT)),category(row.get(AD_CAT))
            if c1 and c2:cats.append((c1,c2))
        dif=np.array([math.log2(a/b) for a,b in pairs])
        n=len(pairs); hit=int(sum(abs(dif)<=1.00000001)); same=sum(x==y for x,y in cats)
        from scipy.stats import wilcoxon
        pval=float(wilcoxon(dif).pvalue) if n>=2 and np.any(dif!=0) else None
        stats.append([drug,n,hit,100*hit/n if n else None,len(cats),same,100*same/len(cats) if cats else None,float(np.mean(dif)) if n else None,pval])
        if n:plotted.append((drug,pairs))
    if plotted:
        fig,axes=plt.subplots(len(plotted),1,figsize=(6,3.5*len(plotted)),squeeze=False,layout="constrained")
        for ax,(d,pairs) in zip(axes.flat,plotted):
            x,y=np.log2(np.array(pairs)).T
            ax.scatter((x+y)/2,x-y,alpha=.65)
            ax.axhline(np.mean(x-y),color="#076b66")
            if len(x)>1:
                for sign in [-1,1]:ax.axhline(np.mean(x-y)+sign*1.96*np.std(x-y,ddof=1),linestyle="--",color="gray")
            ax.set(title=f"{d} • {len(x)} exact pairs",xlabel="Mean log₂ MIC",ylabel="log₂ rounded E-strip − log₂ AD")
        fig.savefig(figdir/"Paired_EST_AD.png",dpi=150)
        plt.close(fig)
    else:audit.add("paired MIC figure","skipped","no exact clinical E-strip / AD pairs")
    agreement_tables(wb,data)
    for method in (ROUNDED,AD_MIC):method_profiles(data,figdir,method,audit)
    style_workbook(wb)
    wb.save(root/"NG_AMR_MIC.xlsx")

def report_gallery(root):
    for method, sub in (("EST", "EST_rounded"), ("AD", "AD")):
        figures = sorted((root / "MIC_Figures" / sub).glob("*.png"))
        profiles = "EST_rounded_genotype_phenotype" if method == "EST" else "AD_genotype_phenotype"
        figures.extend(sorted((root / profiles).glob("*.png")))
        if figures:
            target = root / "MIC_Figures" / sub
            target.mkdir(exist_ok=True)
            if not build_mic_illustration_pptx(str(target), [str(p) for p in figures], str(root)):
                raise RuntimeError(f"Could not create {method} presentation")
