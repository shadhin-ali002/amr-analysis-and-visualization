"""MIC workbook: isolate MICs, method details, QC, paired and categorical agreement."""
from __future__ import annotations

import math
import re

import numpy as np
import openpyxl

from ng_amr.core.constants import (
    AD_CAT,
    AD_MCF,
    AD_MIC,
    ANTIBIOTIC_ORDER,
    CATEGORIES,
    DATE_FIELDS,
    EST_CAT,
    EST_MCF,
    MIC_FIELDS,
    QC,
    ROUNDED,
    WGS_CAT,
)
from ng_amr.io.workbook import (
    measurement,
)


def style_workbook(wb):
    from openpyxl.styles import Alignment, Font, PatternFill
    for sheet in wb:
        sheet.freeze_panes="B2"
        sheet.auto_filter.ref=sheet.dimensions
        sheet.row_dimensions[1].height=38
        for cell in sheet[1]:
            cell.font=Font(color="FFFFFF",bold=True)
            cell.fill=PatternFill("solid",fgColor="234D62")
            cell.alignment=Alignment(wrap_text=True,vertical="center")
            sheet.column_dimensions[cell.column_letter].width=min(52,max(17,len(str(cell.value))*.75))


def category(value):
    return {"s":"S","sensitive":"S","susceptible":"S","i":"I","intermediate":"I","r":"R","resistant":"R"}.get(str(value).strip().lower())


def agreement_tables(wb,book):
    """Recorded categories with explicit method pairs and denominators."""
    from itertools import combinations

    from scipy.stats import binomtest
    ws=wb.create_sheet("Method pairs")
    ws.append(["Antibiotic","Method 1","Method 2","Paired S/I/R","Same S/I/R","Agreement (%)",
               "Paired R or NS binary","Only method 1 R or NS","Only method 2 R or NS","Exact McNemar p"])
    ref=wb.create_sheet("AD comparator errors")
    ref.append(["Antibiotic","S/I/R pairs","Major errors (EST R, AD S)","AD S denominator","Major errors (%)",
                "Very major errors (EST S, AD R)","AD R denominator","Very major errors (%)","Minor errors (I versus S or R)","Minor errors (%)"])
    for drug in ANTIBIOTIC_ORDER:
        records=[r for (k,d),r in book.values.items() if d==drug and k!=QC]
        for left,right in combinations(CATEGORIES,2):
            binary=[];pairs=[]
            for row in records:
                a,b=row.get(left),row.get(right)
                x,y=category(a),category(b)
                if x and y:pairs.append((x,y))
                def resistant(value):
                    k=re.sub(r"[^a-z]","",str(value).lower())
                    return True if k in {"ns","nonsusceptible"} else None if not category(value) else category(value)=="R"
                x,y=resistant(a),resistant(b)
                if x is not None and y is not None:binary.append((x,y))
            n=len(pairs);matches=sum(x==y for x,y in pairs)
            a=sum(x and not y for x,y in binary);b=sum(y and not x for x,y in binary)
            ws.append([drug,left,right,n,matches,100*matches/n if n else None,len(binary),a,b,
                       binomtest(a,a+b,.5).pvalue if a+b else 1 if binary else None])
            if left==EST_CAT and right==AD_CAT:
                major=sum(x=="R" and y=="S" for x,y in pairs);vmajor=sum(x=="S" and y=="R" for x,y in pairs)
                sd=sum(y=="S" for _,y in pairs);rd=sum(y=="R" for _,y in pairs)
                minor=sum((x=="I")!=(y=="I") for x,y in pairs)
                ref.append([drug,n,major,sd,100*major/sd if sd else None,vmajor,rd,100*vmajor/rd if rd else None,minor,100*minor/n if n else None])


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
    stats=wb.create_sheet("Paired agreement")
    stats.append(["Antibiotic","Exact quantitative pairs","Within one doubling dilution","Essential agreement (%)","Category pairs","Identical categories","Categorical agreement (%)","Log2 bias EST minus AD","Wilcoxon paired p"])
    concord=wb.create_sheet("Recorded WGS agreement")
    concord.append(["Antibiotic","Phenotype method","Paired clinical categories","Concordant","Agreement (%)"])
    for method in (ROUNDED, AD_MIC):
        cat=EST_CAT if method==ROUNDED else AD_CAT
        if not any(key!=QC and row.get(method) is not None for (key,d),row in data.values.items()):
            audit.add(method,"skipped","no clinical MICs in the source field")
            continue
        for d in ANTIBIOTIC_ORDER:
            pairs=[(category(row.get(cat)),category(row.get(WGS_CAT))) for (k,drug),row in data.values.items()
                   if k!=QC and drug==d and category(row.get(cat)) and category(row.get(WGS_CAT))]
            agree=sum(a==b for a,b in pairs)
            concord.append([d,method,len(pairs),agree,100*agree/len(pairs) if pairs else None])
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
    agreement_tables(wb,data)
    style_workbook(wb)
    wb.save(root/"NG_AMR_MIC.xlsx")
