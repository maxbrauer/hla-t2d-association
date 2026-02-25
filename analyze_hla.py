#!/usr/bin/env python3
"""
HLA-T2D Association Analyzer
==============================
Analyzes HLA allele data for Type 2 Diabetes risk assessment.
Supports classical HLA typing input and SNP-level analysis.

Usage:
    python analyze_hla.py --demo
    python analyze_hla.py --input hla_typing.csv
    python analyze_hla.py --alleles "DRB1*03:01,DRB1*04:01,DQB1*03:02,DQB1*02:01"
"""

import json
import csv
import argparse
import sys
from pathlib import Path
from dataclasses import dataclass, field, asdict
from collections import defaultdict

DB_PATH = Path(__file__).parent / "data" / "hla_t2d_db.json"


def load_database(db_path: Path = DB_PATH) -> dict:
    with open(db_path, "r", encoding="utf-8") as f:
        return json.load(f)


@dataclass
class HLATyping:
    alleles: list  # List of allele strings e.g. ["DRB1*03:01", "DRB1*15:01", ...]
    snps: dict = field(default_factory=dict)  # {rsid: genotype}


@dataclass
class AlleleResult:
    allele: str
    gene: str
    hla_class: str
    association: str
    odds_ratio: float
    p_value: float
    description: str
    mechanism: str
    haplotype: str
    evidence_level: str
    copies: int  # 0, 1, or 2


@dataclass
class HaplotypeResult:
    haplotype: str
    risk_category: str
    odds_ratio: float
    phenotype: str
    description: str
    detected: bool


@dataclass
class RiskReport:
    overall_risk: str  # "very_high", "high", "moderate", "low", "protective"
    combined_odds_ratio: float
    risk_alleles: list
    protective_alleles: list
    neutral_alleles: list
    haplotype_results: list
    snp_results: list
    alleles_analyzed: int
    alleles_found_in_db: int
    mechanisms_involved: list
    recommendations: list


def parse_alleles_string(alleles_str: str) -> HLATyping:
    """Parse comma-separated HLA alleles."""
    alleles = [a.strip().upper() for a in alleles_str.split(",") if a.strip()]
    # Normalize: ensure HLA- prefix format consistency
    normalized = []
    for a in alleles:
        a = a.replace("HLA-", "")
        normalized.append(a)
    return HLATyping(alleles=normalized)


def parse_csv_input(filepath: str) -> HLATyping:
    """Parse CSV with columns: locus, allele1, allele2 (or allele column)."""
    alleles = []
    snps = {}
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row_lower = {k.lower().strip(): v.strip() for k, v in row.items()}
            # HLA allele rows
            if "allele" in row_lower:
                a = row_lower["allele"].upper().replace("HLA-", "")
                if a:
                    alleles.append(a)
            elif "allele1" in row_lower:
                for key in ["allele1", "allele2"]:
                    a = row_lower.get(key, "").upper().replace("HLA-", "")
                    if a:
                        alleles.append(a)
            # SNP rows
            if "rsid" in row_lower and "genotype" in row_lower:
                snps[row_lower["rsid"]] = row_lower["genotype"]
    return HLATyping(alleles=alleles, snps=snps)


def generate_demo_typing() -> HLATyping:
    """Generate a demo HLA typing for testing."""
    return HLATyping(
        alleles=[
            "DRB1*03:01", "DRB1*15:01",
            "DQA1*05:01", "DQA1*01:02",
            "DQB1*02:01", "DQB1*06:02",
            "A*02:01", "A*24:02",
            "B*39:06", "B*07:02",
            "C*01:02", "C*07:02",
            "DPB1*04:01", "DPB1*02:01",
        ],
        snps={
            "rs9272346": "GA",
            "rs2647044": "CT",
            "rs3129889": "AG",
            "rs9275596": "TC",
            "rs1063355": "CT",
            "rs2395185": "TG",
        }
    )


def analyze_hla(typing: HLATyping, db: dict) -> RiskReport:
    """Run full HLA-T2D association analysis."""
    risk_alleles = []
    protective_alleles = []
    neutral_alleles = []
    found = 0
    combined_log_or = 0.0

    # Match alleles against database
    for allele_entry in db["allele_associations"]:
        db_allele = allele_entry["allele"].replace("HLA-", "")
        copies = sum(1 for a in typing.alleles if a == db_allele)
        if copies == 0:
            continue
        found += 1
        result = AlleleResult(
            allele=allele_entry["allele"],
            gene=allele_entry["gene"],
            hla_class=allele_entry["class"],
            association=allele_entry["association"],
            odds_ratio=allele_entry["odds_ratio"],
            p_value=allele_entry["p_value"],
            description=allele_entry["description"],
            mechanism=allele_entry["mechanism"],
            haplotype=allele_entry["haplotype"],
            evidence_level=allele_entry["evidence_level"],
            copies=copies,
        )
        # Accumulate log(OR) weighted by copies
        import math
        combined_log_or += math.log(allele_entry["odds_ratio"]) * copies

        if allele_entry["association"] == "risk":
            risk_alleles.append(asdict(result))
        elif allele_entry["association"] == "protective":
            protective_alleles.append(asdict(result))
        else:
            neutral_alleles.append(asdict(result))

    # Haplotype detection
    haplotype_results = []
    allele_set = set(typing.alleles)
    for hap in db["haplotype_risks"]:
        # Check if haplotype components are present
        hap_name = hap["haplotype"]
        detected = False
        # Simple detection: check if key alleles from haplotype name are in typing
        if "DRB1*03:01" in allele_set and "DRB1*04:01" in allele_set and "DR3/DR4" in hap_name:
            detected = True
        elif "DRB1*04:01" in allele_set and "DQB1*03:02" in allele_set and "DR4-DQ8" in hap_name:
            detected = True
        elif "DRB1*03:01" in allele_set and "DQB1*02:01" in allele_set and "DR3-DQ2" in hap_name:
            detected = True
        elif "DRB1*09:01" in allele_set and "DQB1*03:03" in allele_set and "DR9-DQ9" in hap_name:
            detected = True
        elif "DRB1*15:01" in allele_set and "DQB1*06:02" in allele_set and "DR15-DQ6" in hap_name:
            detected = True
        elif "DRB1*07:01" in allele_set and "DQB1*02:02" in allele_set and "DR7-DQ2" in hap_name:
            detected = True
        elif "DRB1*13:01" in allele_set and "DQB1*06:03" in allele_set and "DR13-DQ6" in hap_name:
            detected = True

        haplotype_results.append(asdict(HaplotypeResult(
            haplotype=hap["haplotype"],
            risk_category=hap["risk_category"],
            odds_ratio=hap["odds_ratio"],
            phenotype=hap["phenotype"],
            description=hap["description"],
            detected=detected,
        )))

    # SNP analysis
    snp_results = []
    for snp in db["snp_associations"]:
        if snp["rsid"] in typing.snps:
            geno = typing.snps[snp["rsid"]]
            risk_count = geno.count(snp["alt"])
            snp_results.append({
                **snp,
                "genotype": geno,
                "risk_allele_count": risk_count,
            })
            import math
            combined_log_or += math.log(snp["odds_ratio"]) * risk_count

    # Overall risk assessment
    import math
    combined_or = math.exp(combined_log_or) if combined_log_or != 0 else 1.0

    detected_risk_haplotypes = [h for h in haplotype_results if h["detected"] and h["risk_category"] in ("very_high", "high")]
    detected_protective = [h for h in haplotype_results if h["detected"] and "protective" in h["risk_category"]]

    if detected_risk_haplotypes and not detected_protective:
        overall = "very_high" if any(h["risk_category"] == "very_high" for h in detected_risk_haplotypes) else "high"
    elif combined_or > 2.0:
        overall = "high"
    elif combined_or > 1.3:
        overall = "moderate"
    elif combined_or < 0.6:
        overall = "protective"
    elif combined_or < 0.85:
        overall = "low"
    else:
        overall = "moderate"

    if detected_protective and not detected_risk_haplotypes:
        overall = "protective"

    # Mechanisms involved
    mechanisms = []
    if any(a["hla_class"] == "II" for a in risk_alleles):
        mechanisms.append("antigen_presentation")
    if any("DQB1" in a["allele"] for a in risk_alleles + protective_alleles):
        mechanisms.append("asp57_effect")
    if len(risk_alleles) >= 2 and any("DQA1" in al for al in typing.alleles):
        mechanisms.append("trans_complementation")
    if any(a["hla_class"] == "I" for a in risk_alleles):
        mechanisms.append("nk_cell_regulation")

    mechanism_details = [m for m in db["biological_mechanisms"] if m["id"] in mechanisms]

    # Recommendations
    recommendations = []
    if overall in ("very_high", "high"):
        recommendations.extend([
            "Consider GAD65 and IA-2 autoantibody testing to rule out LADA",
            "Monitor C-peptide levels to assess beta-cell function trajectory",
            "HbA1c monitoring every 3 months during first year of diagnosis",
            "Early insulin consideration if rapid beta-cell decline detected",
            "Screen first-degree relatives for autoimmune diabetes markers",
        ])
    elif overall == "moderate":
        recommendations.extend([
            "Standard T2D management with attention to autoimmune markers",
            "Consider GAD antibody testing if atypical presentation",
            "Lifestyle intervention: Mediterranean diet, regular exercise",
            "Regular metabolic panel monitoring",
        ])
    else:
        recommendations.extend([
            "Favorable HLA profile for T2D risk",
            "Standard preventive measures: healthy weight, regular exercise",
            "Routine metabolic screening per age-appropriate guidelines",
        ])

    return RiskReport(
        overall_risk=overall,
        combined_odds_ratio=round(combined_or, 3),
        risk_alleles=risk_alleles,
        protective_alleles=protective_alleles,
        neutral_alleles=neutral_alleles,
        haplotype_results=haplotype_results,
        snp_results=snp_results,
        alleles_analyzed=len(typing.alleles),
        alleles_found_in_db=found,
        mechanisms_involved=mechanism_details,
        recommendations=recommendations,
    )


def print_report(report: RiskReport):
    """Print formatted report to console."""
    risk_labels = {
        "very_high": "VERY HIGH RISK",
        "high": "HIGH RISK",
        "moderate": "MODERATE RISK",
        "low": "LOW RISK",
        "protective": "PROTECTIVE PROFILE",
    }

    print("\n" + "=" * 70)
    print("     HLA-T2D ASSOCIATION ANALYSIS REPORT")
    print("=" * 70)
    print(f"\n  Alleles analyzed: {report.alleles_analyzed}")
    print(f"  Alleles in database: {report.alleles_found_in_db}")
    print(f"  Combined Odds Ratio: {report.combined_odds_ratio:.3f}")
    print(f"\n  OVERALL RISK: {risk_labels.get(report.overall_risk, report.overall_risk)}")

    if report.risk_alleles:
        print("\n" + "-" * 70)
        print("  RISK ALLELES")
        print("-" * 70)
        for a in report.risk_alleles:
            print(f"  {a['allele']:20s} OR={a['odds_ratio']:.2f}  x{a['copies']} copy  [{a['evidence_level']}]")
            print(f"    {a['description'][:90]}...")

    if report.protective_alleles:
        print("\n" + "-" * 70)
        print("  PROTECTIVE ALLELES")
        print("-" * 70)
        for a in report.protective_alleles:
            print(f"  {a['allele']:20s} OR={a['odds_ratio']:.2f}  x{a['copies']} copy  [{a['evidence_level']}]")
            print(f"    {a['description'][:90]}...")

    detected_haps = [h for h in report.haplotype_results if h["detected"]]
    if detected_haps:
        print("\n" + "-" * 70)
        print("  DETECTED HAPLOTYPES")
        print("-" * 70)
        for h in detected_haps:
            print(f"  {h['haplotype']}")
            print(f"    Risk: {h['risk_category'].upper()} | OR={h['odds_ratio']:.2f} | {h['phenotype']}")
            print(f"    {h['description'][:90]}...")

    if report.recommendations:
        print("\n" + "-" * 70)
        print("  RECOMMENDATIONS")
        print("-" * 70)
        for r in report.recommendations:
            print(f"  - {r}")

    print("\n" + "=" * 70)
    print("  Note: For research/educational purposes only.")
    print("  Consult an immunogeneticist for clinical interpretation.")
    print("=" * 70 + "\n")


def main():
    parser = argparse.ArgumentParser(description="HLA-T2D Association Analyzer")
    parser.add_argument("--input", "-i", help="Input CSV file with HLA typing")
    parser.add_argument("--alleles", "-a", help="Comma-separated HLA alleles")
    parser.add_argument("--demo", action="store_true", help="Run with demo data")
    parser.add_argument("--output", "-o", help="Output JSON path")
    parser.add_argument("--db", default=str(DB_PATH), help="Database path")

    args = parser.parse_args()

    if not args.demo and not args.input and not args.alleles:
        parser.print_help()
        print("\nError: Provide --alleles, --input, or --demo")
        sys.exit(1)

    db = load_database(Path(args.db))
    print(f"Loaded {len(db['allele_associations'])} HLA allele associations")
    print(f"Loaded {len(db['snp_associations'])} SNP associations")
    print(f"Loaded {len(db['haplotype_risks'])} haplotype risk profiles")

    if args.demo:
        typing = generate_demo_typing()
        print("Running DEMO analysis...")
    elif args.alleles:
        typing = parse_alleles_string(args.alleles)
    else:
        typing = parse_csv_input(args.input)

    print(f"Input: {len(typing.alleles)} alleles, {len(typing.snps)} SNPs")

    report = analyze_hla(typing, db)
    print_report(report)

    output = args.output or str(Path(__file__).parent / "data" / "hla_report.json")
    with open(output, "w", encoding="utf-8") as f:
        json.dump(asdict(report), f, indent=2, ensure_ascii=False)
    print(f"Report exported to {output}")


if __name__ == "__main__":
    main()
