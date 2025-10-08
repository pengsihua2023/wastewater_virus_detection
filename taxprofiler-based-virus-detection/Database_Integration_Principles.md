# Fundamental Principles of Kraken2 Dual-Database Integration

## 📖 Core Concept

Using two different databases (e.g., RVDB and NCBI Viral RefSeq) for Kraken2 classification and then integrating the results aims to **balance sensitivity and specificity** to obtain more reliable and comprehensive viral detection results.

## 🎯 Five Fundamental Principles

### Principle 1: Confidence Stratification Principle

**Core Idea**: Classify results into different confidence levels based on the strength of detection evidence

#### Stratification Framework

```
┌─────────────────────────────────────────────────────────┐
│                  Confidence Pyramid                     │
├─────────────────────────────────────────────────────────┤
│                                                         │
│            High Confidence                              │
│         ✓ Intersection + High reads                     │
│              ↓ Direct reporting                         │
│                                                         │
│     ┌───────────────────────────────────┐               │
│     │   Medium Confidence               │               │
│     │  ! Single DB + High reads         │               │
│     │      ↓ Validation required        │               │
│     └───────────────────────────────────┘               │
│                                                         │
│  ┌──────────────────────────────────────────┐           │
│  │    Low Confidence                        │           │
│  │   ✗ Single DB + Low reads                │           │
│  │       ↓ Potential false positive         │           │
│  └──────────────────────────────────────────┘           │
└─────────────────────────────────────────────────────────┘
```

#### Implementation Details

**High Confidence**:
- Requirement: Detected in both databases
- Read threshold: Direct assigned reads ≥ 100
- Biological meaning: Core virome component
- Reporting strategy: Primary findings, direct reporting

**Medium Confidence**:
- Requirement: Detected in only one database
- Read threshold: Direct assigned reads ≥ 100
- Biological meaning: Potentially rare virus or database-specific detection
- Reporting strategy: Candidate findings, validation required

**Low Confidence**:
- Condition: Single database detection with reads < 50
- Biological meaning: Potential false positive, contamination, or technical noise
- Reporting strategy: Not recommended for reporting unless specifically justified

---

### Principle 2: Direct Read Assignment Principle

**Core Idea**: Count only classifications with direct read assignments to avoid redundant parent-level counting

#### Problem Explanation

Kraken2 reports contain two types of read counts:
- **Total reads**: Sum of reads for this classification and all its descendants
- **Direct reads**: Reads truly assigned to this taxonomic level

#### Incorrect Approach ❌

```
Counting all classifications with total reads:

Kingdom: Viruses        (1000 total reads, 0 direct reads)    ← Counted
  Phylum: Uroviricota   (800 total reads, 0 direct reads)     ← Counted
    Family: Adenoviridae (600 total reads, 0 direct reads)    ← Counted
      Species: Adenovirus C (600 total reads, 600 direct reads) ← Counted

Result: 4 classifications counted, but only 1 actual detection
```

#### Correct Approach ✅

```
Count only classifications with direct reads:

Kingdom: Viruses        (1000 total, 0 direct)    ← Ignored (no direct reads)
  Phylum: Uroviricota   (800 total, 0 direct)     ← Ignored (no direct reads)
    Family: Adenoviridae (600 total, 0 direct)    ← Ignored (no direct reads)
      Species: Adenovirus C (600 total, 600 direct) ← Counted (has direct reads)

Result: 1 classification counted, reflecting true detection
```

#### Implementation Method

```python
# When parsing Kraken2 report
if reads_direct > 0:  # Only retain classifications with direct reads
    results[taxon_name] = {
        'reads_direct': reads_direct,
        'reads_total': reads_total,
        'rank': rank
    }
```

---

### Principle 3: Database Complementarity Principle

**Core Idea**: Leverage different characteristics of two databases for more comprehensive results

#### RVDB Characteristics

**Strengths**:
- 📚 Broad coverage: Contains 5+ million sequences
- 🦠 High diversity: Includes uncultured and environmental viruses
- 🔬 High sensitivity: Suitable for novel virus discovery

**Weaknesses**:
- ⚠️ High redundancy: Multiple sequence versions for same virus
- ⚠️ Variable quality: Contains partially annotated sequences
- ⚠️ Potentially lower specificity: Increased K-mer collisions

#### NCBI Viral RefSeq Characteristics

**Strengths**:
- ✓ High quality: Strict quality control and manual curation
- ✓ Non-redundant: Only representative sequences per species
- ✓ Standardized annotation: Accurate taxonomic information
- ✓ High specificity: Low false positive rate

**Weaknesses**:
- ✗ Limited coverage: Only known, validated viruses
- ✗ Relatively lower sensitivity: Weaker novel virus detection

#### Complementarity Strategy

```
Decision Tree by Sample Type:

Clinical Sample (Known pathogen detection)
    ↓
Prioritize NCBI RefSeq results
    ├─ High confidence: Intersection
    └─ Medium confidence: NCBI unique

Environmental Sample (Viral diversity study)
    ↓
Fully utilize RVDB sensitivity
    ├─ High confidence: Intersection
    └─ Medium confidence: RVDB unique (potential novel viruses)

Comprehensive Study (Discovery + Identification)
    ↓
Combine strengths of both
    ├─ High confidence: Intersection (core virome)
    ├─ RVDB unique: Potential novel viruses (validation needed)
    └─ NCBI unique: Accurately identified known viruses
```

---

### Principle 4: Threshold Setting Principle

**Core Idea**: Set read thresholds reasonably based on research objectives and sample types

#### Threshold Types

**1. High Confidence Threshold**

Default value: 100 direct reads

Adjustment criteria:
- High sequencing depth (>10M reads) → Increase to 150-200
- Low sequencing depth (<1M reads) → Decrease to 50-75
- Clinical diagnostics → Can decrease (more sensitive)
- Environmental samples → Can increase (reduce noise)

**2. Medium Confidence Threshold**

Default value: 50 direct reads

Adjustment criteria:
- Discovery-oriented research → Decrease to 25-30
- Conservative analysis → Increase to 75-100

**3. Minimum Direct Reads Threshold**

Default value: 0 (retain all classifications with direct reads)

Adjustment criteria:
- Excessive noise → Set to 5-10
- Maximum sensitivity → Keep at 0

#### Threshold Setting Examples

```python
# Clinical sample (pathogen diagnosis)
reads_threshold_high = 50      # More sensitive
reads_threshold_medium = 25
min_direct_reads = 5

# Standard research sample
reads_threshold_high = 100     # Balanced
reads_threshold_medium = 50
min_direct_reads = 0

# Environmental sample (high diversity)
reads_threshold_high = 200     # More stringent
reads_threshold_medium = 100
min_direct_reads = 10

# Deep sequencing sample (>50M reads)
reads_threshold_high = 500     # Correspondingly higher
reads_threshold_medium = 200
min_direct_reads = 20
```

#### Dynamic Threshold Strategy

```
Automatically adjust based on total reads:

if total_reads < 1M:
    high_threshold = 50
elif total_reads < 10M:
    high_threshold = 100
elif total_reads < 50M:
    high_threshold = 200
else:
    high_threshold = 500
```

---

### Principle 5: Validation and Reporting Principle

**Core Idea**: Apply different validation and reporting strategies for different confidence levels

#### High Confidence Results

**Validation requirement**: Minimal

**Reporting strategy**:
- ✓ Report directly in main text
- ✓ Present as primary findings
- ✓ Suitable for biological interpretation
- ✓ Use for quantitative analysis

**Example phrasing**:
```
"We identified 15 high-confidence viral species, including 
Human betaherpesvirus 5 (CMV, 12.5% of viral reads), 
Epstein-Barr virus (5.2%), and Torque teno virus (3.1%)."
```

#### Medium Confidence Results

**Validation requirement**: Strongly recommended

**Validation methods**:
1. **Sequence validation**
   ```bash
   # Extract reads matching the virus
   # BLAST validation of reads
   blastn -query candidate_virus_reads.fasta -db nt -max_target_seqs 5
   ```

2. **Coverage analysis**
   ```bash
   # Check viral genome coverage
   # Coverage >1% and uniform distribution → More reliable
   ```

3. **Sequence assembly**
   ```bash
   # If sufficient reads (>1000), attempt assembly
   # Validate assembled contigs
   ```

**Reporting strategy**:
- Briefly mention in results section
- Details in supplementary materials
- Label as "candidate" or "pending validation"

**Example phrasing**:
```
"Additionally, 8 candidate viruses were detected exclusively in 
RVDB with >100 reads (Supplementary Table S2), potentially 
representing rare or recently discovered viruses requiring 
further validation."
```

#### Low Confidence Results

**Validation requirement**: Optional

**Reporting strategy**:
- Usually not reported
- Can briefly mention in discussion
- Used to illustrate detection limits or database differences

**Example phrasing**:
```
"A total of 150 low-confidence detections (<50 direct reads) 
were observed, primarily in RVDB, likely representing technical 
noise or database artifacts rather than true viral presence."
```

---

## 🔄 Integration Workflow

### Step 1: Data Preprocessing

```
Kraken2 classification results
    ↓
Parse report files
    ↓
Filter: Only retain classifications with direct_reads > 0
    ↓
Group by taxonomic rank (Species, Genus, Family, etc.)
```

### Step 2: Database Matching

```
RVDB classification list    NCBI classification list
    ↓                              ↓
          Match by taxon name
                ↓
        ┌───────┴───────┐
        │               │
  Intersection      DB-unique
        │               │
        ↓               ↓
  High confidence   Medium/Low confidence
   candidates          candidates
```

### Step 3: Confidence Scoring

```
For each detected taxon:

if (in_RVDB AND in_NCBI):
    if (max_direct_reads >= threshold_high):
        confidence = "High"
    else:
        confidence = "High-Low"
        
else:  # Detected in only one database
    if (max_direct_reads >= threshold_high):
        confidence = "Medium"
    elif (max_direct_reads >= threshold_medium):
        confidence = "Medium-Low"
    else:
        confidence = "Low"
```

### Step 4: Stratified Output

```
High confidence results → Primary findings file
    ↓
Medium confidence results → Candidate viruses file
    ↓
Low confidence results → Complete results file (reference)
    ↓
Statistical summary → Summary report
```

### Step 5: Visualization and Reporting

```
Venn diagram → Show database overlap
    ↓
Bar chart → Top virus abundance
    ↓
Heatmap → Inter-sample comparison
    ↓
Publication tables and figures
```

---

## 📊 Practical Application Cases

### Case 1: Clinical Blood Sample

**Background**: Patient with suspected viral infection, sequencing depth: 5M reads

**Results**:
- RVDB: 35 species
- NCBI: 12 species
- Intersection: 10 species

**Integrated analysis**:
```
High confidence (10, intersection):
- Human cytomegalovirus (CMV): 8.5%
- Epstein-Barr virus (EBV): 3.2%
→ Confirmed infection, direct reporting

Medium confidence (8, single DB high reads):
- Torque teno virus (RVDB unique): 2.1%
- BLAST validation needed → Validation passed → Report as co-infection

Low confidence (29):
- Various low-abundance detections
→ Not reported, likely database noise
```

### Case 2: Seawater Metagenome

**Background**: Environmental viral diversity study, sequencing depth: 50M reads

**Results**:
- RVDB: 280 species
- NCBI: 45 species
- Intersection: 38 species

**Integrated analysis**:
```
High confidence (38, intersection):
- Core marine viral community
- Includes known phages and algal viruses
→ Primary findings, detailed reporting

Medium confidence (RVDB unique 242):
- 150 are "uncultured marine viruses"
- Represent true environmental viral diversity
→ Report total count, detailed description after partial validation

Medium confidence (NCBI unique 7):
- May be unclassified in RVDB due to K-mer conflicts
- Should be true viruses
→ Add to core virome reporting
```

---

## 🎯 Key Decision Points

### Decision 1: Use One or Two Databases?

```
Known pathogen detection only → NCBI RefSeq sufficient

Virus discovery research → Must use RVDB

Comprehensive study → Combine both (recommended)
```

### Decision 2: How to Set Thresholds?

```
Based on data:
- Deep sequencing → Increase thresholds
- Shallow sequencing → Decrease thresholds

Based on objectives:
- Discovery-oriented → Lower thresholds (more sensitive)
- Diagnostic-oriented → Balanced thresholds
- Quantitative study → Higher thresholds (more accurate)

Based on samples:
- Clinical samples → Can lower thresholds
- Environmental samples → Can increase thresholds (more noise)
```

### Decision 3: How to Report Results?

```
High confidence (intersection):
→ Detailed reporting in main text

Medium confidence (single DB high reads):
→ Brief reporting + supplementary materials

Low confidence (single DB low reads):
→ Usually not reported, or mentioned only in methods
```

---

## 💡 Best Practice Recommendations

### 1. Transparency Principle
- ✓ Clearly state database versions and download dates
- ✓ Report threshold settings and rationale
- ✓ Distinguish results of different confidence levels

### 2. Reproducibility Principle
- ✓ Save complete configuration files
- ✓ Record all parameters
- ✓ Provide original result files

### 3. Conservative Principle
- ✓ Be cautious with single-database detections
- ✓ Independently validate important findings
- ✓ Choose more conservative interpretation when uncertain

### 4. Biological Plausibility Principle
- ✓ Consider sample type and source
- ✓ Reference known epidemiological data
- ✓ Check taxonomic reasonableness

### 5. Statistical Rigor Principle
- ✓ Consider false positive rates
- ✓ Multiple testing correction (if needed)
- ✓ Report confidence intervals

---

## 📚 Theoretical Foundation

### Bayesian Perspective

The intersection of two databases can be viewed as an improvement in **Bayesian posterior probability**:

```
P(Virus present | Detected in both DBs) > P(Virus present | Detected in one DB)

Prior: Database quality and coverage
Likelihood: Read quantity and quality
Posterior: Integrated confidence
```

### Set Theory Perspective

```
Universe U: All possible viruses
RVDB set R: Viruses detected by RVDB
NCBI set N: Viruses detected by NCBI

High confidence = R ∩ N (intersection)
Medium confidence = (R ∪ N) - (R ∩ N) (union minus intersection, i.e., symmetric difference)
True virome ≈ R ∩ N + subset of (R △ N)
```

### Sensitivity-Specificity Balance

```
RVDB: High sensitivity, moderate specificity
NCBI: Moderate sensitivity, high specificity
Intersection: Moderate sensitivity, high specificity (balance point)
```

---

## 🔚 Summary

The core of integrating two Kraken2 database results is:

1. **Stratified management**: Assign confidence based on evidence strength
2. **Accurate counting**: Count only direct reads to avoid duplication
3. **Complementary utilization**: Combine advantages of both databases
4. **Reasonable thresholds**: Adjust standards based on specific situations
5. **Rigorous reporting**: Distinguish findings at different levels

By following these principles, you can obtain **comprehensive yet reliable** viral classification results, maximizing true positive detection rate while minimizing false positive risk.

---

**Document Version**: 1.0  
**Creation Date**: October 7, 2025  
**Applicable Tools**: integrate_66ce4dde_EN.py / integrate_66ce4dde.py




