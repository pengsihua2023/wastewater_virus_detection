#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kraken2 Results Integration Script - For sample 66ce4dde (Improved v2)
Only retains classifications with direct read assignments, avoiding parent-level duplication

Input files: 
  - 66ce4dde_kraken2_RVDB.txt
  - 66ce4dde_kraken2_NCBI.txt

Author: Cursor AI Assistant
Date: 2025-10-07
Version: 2.0
"""

import pandas as pd
import sys
import os
from collections import defaultdict

def parse_kraken_report(report_file, min_direct_reads=0):
    """
    Parse Kraken2 report file, retaining only classifications with direct read assignments
    
    Parameters:
        report_file: Path to Kraken2 report file
        min_direct_reads: Minimum direct reads (default 0, retains all with direct reads)
    
    Returns:
        Dictionary containing all classifications with direct read assignments
    """
    print(f"Parsing: {report_file}")
    results = {}
    
    try:
        with open(report_file, 'r', encoding='utf-8') as f:
            for line in f:
                # Kraken2 report format:
                # Percentage  Total_reads  Direct_reads  Rank  NCBI_TaxID  Name
                fields = line.strip().split('\t')
                
                if len(fields) >= 6:
                    percent = float(fields[0].strip())
                    reads_total = int(fields[1].strip())
                    reads_direct = int(fields[2].strip())
                    rank = fields[3].strip()
                    taxid = fields[4].strip()
                    name = fields[5].strip()
                    
                    # Only retain classifications with direct read assignments
                    # This avoids counting parent-level summaries (Kingdom, Phylum, etc.)
                    if reads_direct > min_direct_reads:
                        results[name] = {
                            'percent': percent,
                            'reads_total': reads_total,
                            'reads_direct': reads_direct,
                            'rank': rank,
                            'taxid': taxid
                        }
        
        print(f"  - Parsed {len(results)} classification entries with direct reads")
        
        # Count by taxonomic rank
        rank_counts = {}
        for data in results.values():
            rank = data['rank']
            rank_counts[rank] = rank_counts.get(rank, 0) + 1
        
        print(f"  - Rank distribution: ", end="")
        rank_names = {'S': 'Species', 'G': 'Genus', 'F': 'Family', 'O': 'Order', 
                     'C': 'Class', 'P': 'Phylum', 'K': 'Kingdom', 'D': 'Domain'}
        for rank in ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'D']:
            if rank in rank_counts:
                print(f"{rank_names.get(rank, rank)}:{rank_counts[rank]} ", end="")
        print()
        
        return results
        
    except Exception as e:
        print(f"Error: Cannot parse file {report_file}")
        print(f"Error message: {e}")
        sys.exit(1)

def get_rank_priority(rank):
    """
    Get taxonomic rank priority (for sorting)
    
    Parameters:
        rank: Rank code
    
    Returns:
        Priority number (lower is more important)
    """
    rank_order = {
        'D': 1,   # Domain
        'K': 2,   # Kingdom  
        'P': 3,   # Phylum
        'C': 4,   # Class
        'O': 5,   # Order
        'F': 6,   # Family
        'G': 7,   # Genus
        'S': 8,   # Species
        'S1': 9,  # Subspecies
        'U': 10   # Unclassified
    }
    return rank_order.get(rank, 99)

def get_rank_name(rank):
    """
    Get full taxonomic rank name
    """
    rank_names = {
        'D': 'Domain',
        'K': 'Kingdom',
        'P': 'Phylum',
        'C': 'Class',
        'O': 'Order',
        'F': 'Family',
        'G': 'Genus',
        'S': 'Species',
        'S1': 'Subspecies',
        'U': 'Unclassified'
    }
    return rank_names.get(rank, rank)

def integrate_results(rvdb_file, ncbi_file, 
                     reads_threshold_high=100, 
                     reads_threshold_medium=50,
                     min_direct_reads=0,
                     output_prefix='66ce4dde_integrated'):
    """
    Integrate classification results from two databases (improved version - only counts direct reads)
    
    Parameters:
        rvdb_file: RVDB database result file
        ncbi_file: NCBI RefSeq database result file
        reads_threshold_high: High confidence reads threshold
        reads_threshold_medium: Medium confidence reads threshold
        min_direct_reads: Minimum direct reads
        output_prefix: Output file prefix
    """
    
    print("\n" + "="*70)
    print("Kraken2 Results Integration Analysis (Improved v2)")
    print("Sample: 66ce4dde")
    print("Improvement: Only counts classifications with direct reads, avoiding parent duplication")
    print("="*70 + "\n")
    
    # Parse both report files (only retain classifications with direct reads)
    rvdb_results = parse_kraken_report(rvdb_file, min_direct_reads)
    ncbi_results = parse_kraken_report(ncbi_file, min_direct_reads)
    
    # Get all detected taxon names
    all_taxa = set(rvdb_results.keys()) | set(ncbi_results.keys())
    print(f"\nTotal detected: {len(all_taxa)} classification entries with direct reads")
    
    # Integration results
    integrated = []
    
    for taxon in all_taxa:
        rvdb_data = rvdb_results.get(taxon, {
            'percent': 0, 
            'reads_total': 0, 
            'reads_direct': 0,
            'rank': 'U',
            'taxid': 'N/A'
        })
        
        ncbi_data = ncbi_results.get(taxon, {
            'percent': 0, 
            'reads_total': 0,
            'reads_direct': 0,
            'rank': 'U',
            'taxid': 'N/A'
        })
        
        # Check if detected in both databases
        in_rvdb = taxon in rvdb_results
        in_ncbi = taxon in ncbi_results
        in_both = in_rvdb and in_ncbi
        
        # Get maximum direct reads (this is the true read count)
        max_direct_reads = max(rvdb_data['reads_direct'], ncbi_data['reads_direct'])
        max_percent = max(rvdb_data['percent'], ncbi_data['percent'])
        
        # Determine confidence level (based on direct reads)
        if in_both:
            if max_direct_reads >= reads_threshold_high:
                confidence = 'High'
                priority = 1
                description = 'Intersection-HighReads'
            else:
                confidence = 'High-Low'
                priority = 2
                description = 'Intersection-LowReads'
        else:
            if max_direct_reads >= reads_threshold_high:
                confidence = 'Medium'
                priority = 3
                source = 'RVDB' if in_rvdb else 'NCBI'
                description = f'{source}Only-HighReads'
            elif max_direct_reads >= reads_threshold_medium:
                confidence = 'Medium-Low'
                priority = 4
                source = 'RVDB' if in_rvdb else 'NCBI'
                description = f'{source}Only-MediumReads'
            else:
                confidence = 'Low'
                priority = 5
                source = 'RVDB' if in_rvdb else 'NCBI'
                description = f'{source}Only-LowReads'
        
        # Get taxonomic rank (prefer the one with data)
        rank = rvdb_data['rank'] if in_rvdb else ncbi_data['rank']
        taxid = rvdb_data['taxid'] if in_rvdb else ncbi_data['taxid']
        
        # Calculate average percentage
        if in_both:
            avg_percent = (rvdb_data['percent'] + ncbi_data['percent']) / 2
        else:
            avg_percent = max_percent
        
        integrated.append({
            'Taxon_Name': taxon,
            'Rank': rank,
            'Rank_Name': get_rank_name(rank),
            'Rank_Priority': get_rank_priority(rank),
            'NCBI_TaxID': taxid,
            'Confidence': confidence,
            'Priority': priority,
            'Description': description,
            'In_Both': 'Yes' if in_both else 'No',
            'In_RVDB': 'Yes' if in_rvdb else 'No',
            'In_NCBI': 'Yes' if in_ncbi else 'No',
            'RVDB_Reads_Direct': rvdb_data['reads_direct'],
            'RVDB_Reads_Total': rvdb_data['reads_total'],
            'RVDB_Percent': rvdb_data['percent'],
            'NCBI_Reads_Direct': ncbi_data['reads_direct'],
            'NCBI_Reads_Total': ncbi_data['reads_total'],
            'NCBI_Percent': ncbi_data['percent'],
            'Max_Direct_Reads': max_direct_reads,
            'Max_Percent': max_percent,
            'Avg_Percent': avg_percent
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(integrated)
    
    # Sort: by priority, then rank, then direct reads
    df = df.sort_values(
        ['Priority', 'Rank_Priority', 'Max_Direct_Reads'], 
        ascending=[True, True, False]
    )
    
    # Save complete results
    full_output = f"{output_prefix}_full.tsv"
    df.to_csv(full_output, sep='\t', index=False, encoding='utf-8')
    print(f"\n✅ Complete results saved to: {full_output}")
    
    # Generate summary statistics
    print("\n" + "="*70)
    print("Integration Summary Statistics (based on direct read classifications)")
    print("="*70)
    
    total_count = len(df)
    high_conf = len(df[df['Confidence'] == 'High'])
    high_low_conf = len(df[df['Confidence'] == 'High-Low'])
    medium_conf = len(df[df['Confidence'] == 'Medium'])
    medium_low_conf = len(df[df['Confidence'] == 'Medium-Low'])
    low_conf = len(df[df['Confidence'] == 'Low'])
    in_both = len(df[df['In_Both'] == 'Yes'])
    
    print(f"\nTotal classifications: {total_count}")
    print(f"  - Detected in both: {in_both} ({in_both/total_count*100:.1f}%)")
    print(f"  - RVDB only: {len(df[df['In_RVDB']=='Yes']) - in_both}")
    print(f"  - NCBI only: {len(df[df['In_NCBI']=='Yes']) - in_both}")
    
    print(f"\nBy confidence level:")
    print(f"  - High confidence: {high_conf}")
    print(f"  - High-Low confidence: {high_low_conf}")
    print(f"  - Medium confidence: {medium_conf}")
    print(f"  - Medium-Low confidence: {medium_low_conf}")
    print(f"  - Low confidence: {low_conf}")
    
    # Statistics by taxonomic rank
    print(f"\nBy taxonomic rank:")
    for rank in ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'D']:
        count = len(df[df['Rank'] == rank])
        if count > 0:
            print(f"  - {get_rank_name(rank)}: {count}")
    
    # Generate high confidence results (recommended for reporting)
    print("\n" + "="*70)
    print("High Confidence Results (Recommended for reporting)")
    print("="*70)
    
    high_confidence = df[df['Confidence'] == 'High'].copy()
    
    if len(high_confidence) > 0:
        # Save species-level high confidence results
        high_conf_species = high_confidence[high_confidence['Rank'] == 'S'].copy()
        
        if len(high_conf_species) > 0:
            # Save high confidence species results
            species_output = f"{output_prefix}_high_confidence_species.tsv"
            high_conf_species.to_csv(species_output, sep='\t', index=False, encoding='utf-8')
            print(f"\n✅ High confidence species results saved to: {species_output}")
            
            # Display top 20
            print(f"\nHigh confidence species (Top 20, sorted by direct reads):")
            print("-" * 130)
            display_cols = ['Taxon_Name', 'RVDB_Reads_Direct', 'NCBI_Reads_Direct', 
                          'RVDB_Percent', 'NCBI_Percent', 'Avg_Percent']
            print(high_conf_species[display_cols].head(20).to_string(index=False))
            
        # Save all-level high confidence results
        all_high_output = f"{output_prefix}_high_confidence_all.tsv"
        high_confidence.to_csv(all_high_output, sep='\t', index=False, encoding='utf-8')
        print(f"\n✅ All-level high confidence results saved to: {all_high_output}")
        
        # Statistics by level
        print(f"\nHigh confidence by taxonomic rank:")
        for rank in ['S', 'G', 'F', 'O', 'C', 'P', 'K']:
            count = len(high_confidence[high_confidence['Rank'] == rank])
            if count > 0:
                print(f"  - {get_rank_name(rank)}: {count}")
    else:
        print("\n⚠️  Warning: No high confidence results detected!")
    
    # Generate candidate virus list (medium confidence)
    print("\n" + "="*70)
    print("Candidate Viruses (Require validation)")
    print("="*70)
    
    candidates = df[df['Confidence'] == 'Medium'].copy()
    
    # Display candidate counts by level
    print(f"\nCandidates by taxonomic rank:")
    for rank in ['S', 'G', 'F', 'O']:
        count = len(candidates[candidates['Rank'] == rank])
        if count > 0:
            print(f"  - {get_rank_name(rank)}: {count}")
    
    candidates_species = candidates[candidates['Rank'] == 'S']
    
    if len(candidates_species) > 0:
        candidates_output = f"{output_prefix}_candidate_viruses.tsv"
        candidates_species.to_csv(candidates_output, sep='\t', index=False, encoding='utf-8')
        print(f"\n✅ Candidate viruses (species-level) saved to: {candidates_output}")
        
        print(f"\nCandidate viruses species (Top 10, sorted by direct reads):")
        print("-" * 130)
        display_cols = ['Taxon_Name', 'In_RVDB', 'In_NCBI', 
                       'Max_Direct_Reads', 'Max_Percent', 'Description']
        print(candidates_species[display_cols].head(10).to_string(index=False))
    else:
        print("\nNo candidate viruses (species-level) detected.")
    
    # Generate Venn diagram data (species-level only)
    print("\n" + "="*70)
    print("Generate Venn Diagram Data (Species-level only)")
    print("="*70)
    
    # Only count species level
    species_df = df[df['Rank'] == 'S']
    rvdb_species = set(species_df[species_df['In_RVDB'] == 'Yes']['Taxon_Name'])
    ncbi_species = set(species_df[species_df['In_NCBI'] == 'Yes']['Taxon_Name'])
    intersection = rvdb_species & ncbi_species
    
    venn_output = f"{output_prefix}_venn_data.txt"
    with open(venn_output, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("Sample 66ce4dde - Venn Diagram Data\n")
        f.write("(Species-level, only classifications with direct reads)\n")
        f.write("="*70 + "\n\n")
        f.write(f"RVDB detected: {len(rvdb_species)} species\n")
        f.write(f"NCBI detected: {len(ncbi_species)} species\n")
        f.write(f"Intersection (high confidence): {len(intersection)} species\n")
        f.write(f"RVDB unique: {len(rvdb_species - ncbi_species)} species\n")
        f.write(f"NCBI unique: {len(ncbi_species - rvdb_species)} species\n\n")
        
        f.write("="*70 + "\n")
        f.write("Intersection Species List (alphabetically sorted)\n")
        f.write("="*70 + "\n")
        for species in sorted(intersection):
            # Add reads information
            species_data = species_df[species_df['Taxon_Name'] == species].iloc[0]
            f.write(f"{species}\n")
            f.write(f"  RVDB: {species_data['RVDB_Reads_Direct']} direct reads ({species_data['RVDB_Percent']:.2f}%)\n")
            f.write(f"  NCBI: {species_data['NCBI_Reads_Direct']} direct reads ({species_data['NCBI_Percent']:.2f}%)\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("RVDB Unique Species List (Top 20, sorted by reads)\n")
        f.write("="*70 + "\n")
        rvdb_only_df = species_df[(species_df['In_RVDB'] == 'Yes') & (species_df['In_NCBI'] == 'No')]
        rvdb_only_df = rvdb_only_df.sort_values('RVDB_Reads_Direct', ascending=False)
        for idx, row in rvdb_only_df.head(20).iterrows():
            f.write(f"{row['Taxon_Name']}\n")
            f.write(f"  RVDB: {row['RVDB_Reads_Direct']} direct reads ({row['RVDB_Percent']:.2f}%)\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("NCBI Unique Species List (Top 20, sorted by reads)\n")
        f.write("="*70 + "\n")
        ncbi_only_df = species_df[(species_df['In_NCBI'] == 'Yes') & (species_df['In_RVDB'] == 'No')]
        ncbi_only_df = ncbi_only_df.sort_values('NCBI_Reads_Direct', ascending=False)
        for idx, row in ncbi_only_df.head(20).iterrows():
            f.write(f"{row['Taxon_Name']}\n")
            f.write(f"  NCBI: {row['NCBI_Reads_Direct']} direct reads ({row['NCBI_Percent']:.2f}%)\n")
    
    print(f"\n✅ Venn diagram data saved to: {venn_output}")
    print(f"\n📊 Species-level statistics (with direct reads):")
    print(f"  - RVDB detected: {len(rvdb_species)} species")
    print(f"  - NCBI detected: {len(ncbi_species)} species")
    print(f"  - Intersection: {len(intersection)} species")
    print(f"  - RVDB unique: {len(rvdb_species - ncbi_species)} species")
    print(f"  - NCBI unique: {len(ncbi_species - rvdb_species)} species")
    
    # Generate summary report
    summary_output = f"{output_prefix}_summary.txt"
    with open(summary_output, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("Sample 66ce4dde Virus Classification Integration Report\n")
        f.write("Improvement: Only counts direct read classifications, avoiding parent duplication\n")
        f.write("="*70 + "\n\n")
        
        f.write("Database Information:\n")
        f.write(f"  - RVDB result file: {os.path.basename(rvdb_file)}\n")
        f.write(f"  - NCBI result file: {os.path.basename(ncbi_file)}\n\n")
        
        f.write("Integration Parameters:\n")
        f.write(f"  - High confidence reads threshold: {reads_threshold_high}\n")
        f.write(f"  - Medium confidence reads threshold: {reads_threshold_medium}\n")
        f.write(f"  - Minimum direct reads: {min_direct_reads}\n\n")
        
        f.write("Integration Statistics (direct read classifications):\n")
        f.write(f"  - Total classifications: {total_count}\n")
        f.write(f"  - Both databases: {in_both}\n")
        f.write(f"  - High confidence: {high_conf}\n")
        f.write(f"  - Medium confidence: {medium_conf}\n")
        f.write(f"  - Low confidence: {low_conf}\n\n")
        
        f.write("Species-level Statistics:\n")
        f.write(f"  - RVDB species count: {len(rvdb_species)}\n")
        f.write(f"  - NCBI species count: {len(ncbi_species)}\n")
        f.write(f"  - Intersection species: {len(intersection)}\n")
        f.write(f"  - RVDB unique: {len(rvdb_species - ncbi_species)}\n")
        f.write(f"  - NCBI unique: {len(ncbi_species - rvdb_species)}\n\n")
        
        if len(high_conf_species) > 0:
            f.write("="*70 + "\n")
            f.write("High Confidence Viruses (Top 10)\n")
            f.write("="*70 + "\n")
            for i, (idx, row) in enumerate(high_conf_species.head(10).iterrows(), 1):
                f.write(f"\n{i}. {row['Taxon_Name']}\n")
                f.write(f"   RVDB: {row['RVDB_Reads_Direct']} direct reads ({row['RVDB_Percent']:.2f}%)\n")
                f.write(f"   NCBI: {row['NCBI_Reads_Direct']} direct reads ({row['NCBI_Percent']:.2f}%)\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("Generated Files List\n")
        f.write("="*70 + "\n")
        f.write(f"1. {full_output} - Complete integration results\n")
        if len(high_confidence) > 0:
            f.write(f"2. {species_output} - High confidence species results\n")
            f.write(f"3. {all_high_output} - All-level high confidence results\n")
        if len(candidates_species) > 0:
            f.write(f"4. {candidates_output} - Candidate virus list\n")
        f.write(f"5. {venn_output} - Venn diagram data\n")
        f.write(f"6. {summary_output} - This summary report\n")
    
    print(f"\n✅ Summary report saved to: {summary_output}")
    
    print("\n" + "="*70)
    print("✅ Integration analysis complete!")
    print("="*70 + "\n")
    
    return df

def main():
    """
    Main function
    """
    # Default file paths (new filenames)
    rvdb_file = "66ce4dde_kraken2_RVDB.txt"
    ncbi_file = "66ce4dde_kraken2_NCBI.txt"
    
    # Use command line arguments if provided
    if len(sys.argv) >= 3:
        rvdb_file = sys.argv[1]
        ncbi_file = sys.argv[2]
    
    # Optional: set thresholds from command line
    reads_threshold_high = 100
    reads_threshold_medium = 50
    min_direct_reads = 0  # Minimum direct reads, 0 = retain all with direct reads
    
    if len(sys.argv) >= 4:
        reads_threshold_high = int(sys.argv[3])
    if len(sys.argv) >= 5:
        reads_threshold_medium = int(sys.argv[4])
    if len(sys.argv) >= 6:
        min_direct_reads = int(sys.argv[5])
    
    # Check if files exist
    if not os.path.exists(rvdb_file):
        print(f"Error: Cannot find RVDB result file: {rvdb_file}")
        print("\nUsage:")
        print("  python integrate_66ce4dde_EN.py [RVDB_file] [NCBI_file] [high_threshold] [medium_threshold] [min_direct_reads]")
        print("\nExample:")
        print("  python integrate_66ce4dde_EN.py \\")
        print("    66ce4dde_kraken2_RVDB.txt \\")
        print("    66ce4dde_kraken2_NCBI.txt \\")
        print("    100 50 0")
        print("\nExplanation:")
        print("  - high_threshold: High confidence reads threshold (default 100)")
        print("  - medium_threshold: Medium confidence reads threshold (default 50)")
        print("  - min_direct_reads: Minimum direct reads (default 0, retains all with direct reads)")
        sys.exit(1)
    
    if not os.path.exists(ncbi_file):
        print(f"Error: Cannot find NCBI result file: {ncbi_file}")
        sys.exit(1)
    
    # Run integration analysis
    integrate_results(rvdb_file, ncbi_file, 
                     reads_threshold_high, 
                     reads_threshold_medium,
                     min_direct_reads)

if __name__ == "__main__":
    main()

