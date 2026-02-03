#!/usr/bin/env python3
"""
Simple script to analyze the last 12 months of malaria simulation data.
- PfPR (Parasite Prevalence Rate)
- Clinical episodes per 1000 population per year

Usage:
    Set db_path variable to your SQLite database file, then run:
    python analyze_last_year.py
"""

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt

def connect_to_db(db_path):
    """Connect to the SQLite database."""
    conn = sqlite3.connect(db_path)
    return conn

def get_last_year_data(conn):
    """Extract data for the last 12 months from the database."""
    query = """
    SELECT 
        md.days_elapsed,
        msd.pfpr_all,
        msd.pfpr_2to10,
        msd.pfpr_under5,
        msd.population,
        msd.clinical_episodes
    FROM monthly_site_data_district msd
    JOIN monthly_data md ON msd.monthly_data_id = md.id
    ORDER BY md.days_elapsed DESC
    LIMIT 12
    """
    df = pd.read_sql_query(query, conn)
    df = df.sort_values('days_elapsed')
    return df

def get_summary_statistics(df):
    """Calculate summary statistics for the last year."""
    # Annual clinical episodes per 1000 = sum(12 months) / mean(population) * 1000
    total_episodes = df['clinical_episodes'].sum()
    mean_population = df['population'].mean()
    annual_episodes_per_1000 = (total_episodes / mean_population) * 1000
    
    return {
        'mean_pfpr_all': df['pfpr_all'].mean(),
        'mean_pfpr2to10': df['pfpr2to10'].mean(),
        'mean_pfpr_under5': df['pfpr_under5'].mean(),
        'mean_population': mean_population,
        'total_clinical_episodes': total_episodes,
        'annual_episodes_per_1000': annual_episodes_per_1000
    }

def print_monthly_data(df):
    """Print monthly data in a formatted table."""
    print("\n" + "="*80)
    print("MONTHLY DATA (Last 12 months)")
    print("="*80)
    print(f"{'Day':<8} {'PfPR All':<10} {'PfPR 2-10':<10} {'PfPR <5':<10} {'Population':<12} {'Episodes':<10}")
    print("-"*80)
    for _, row in df.iterrows():
        print(f"{int(row['days_elapsed']):<8} {row['pfpr_all']:>9.4f} {row['pfpr_2to10']:>9.4f} "
              f"{row['pfpr_under5']:>9.4f} {int(row['population']):>11,} {int(row['clinical_episodes']):>9,}")
    print("="*80)

def print_summary(summary):
    """Print summary statistics."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS (Last 12 months)")
    print("="*80)
    print(f"Mean PfPR (All ages):        {summary['mean_pfpr_all']:>8.2%}")
    print(f"Mean PfPR (2-10 years):      {summary['mean_pfpr2to10']:>8.2%}")
    print(f"Mean PfPR (Under 5):         {summary['mean_pfpr_under5']:>8.2%}")
    print(f"Mean Population:             {summary['mean_population']:>11,.0f}")
    print(f"Total Clinical Episodes:     {summary['total_clinical_episodes']:>11,.0f}")
    print(f"Annual Episodes per 1000:    {summary['annual_episodes_per_1000']:>11,.1f}")
    print("="*80 + "\n")

def plot_data(df):
    """Create visualization of the data."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Panel 1: Clinical episodes per 1000 population
    episodes_per_1000 = (df['clinical_episodes'] / df['population']) * 1000
    axes[0].plot(df['days_elapsed'], episodes_per_1000, marker='o', linewidth=2, markersize=6, color='#e74c3c')
    axes[0].set_xlabel('Days Elapsed')
    axes[0].set_ylabel('Episodes per 1000')
    axes[0].set_title('Clinical Episodes per 1000 Population (Monthly)', fontweight='bold', fontsize=12)
    axes[0].grid(True, alpha=0.3)
    
    # Panel 2: PfPR by age group
    axes[1].plot(df['days_elapsed'], df['pfpr_all'], marker='o', label='All ages', linewidth=2, markersize=6)
    axes[1].plot(df['days_elapsed'], df['pfpr_2to10'], marker='s', label='2-10 years', linewidth=2, markersize=6)
    axes[1].plot(df['days_elapsed'], df['pfpr_under5'], marker='^', label='Under 5', linewidth=2, markersize=6)
    axes[1].set_xlabel('Days Elapsed')
    axes[1].set_ylabel('PfPR')
    axes[1].set_title('Parasite Prevalence Rate by Age Group', fontweight='bold', fontsize=12)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Panel 3: Absolute clinical episodes
    axes[2].plot(df['days_elapsed'], df['clinical_episodes'], marker='o', linewidth=2, markersize=6, color='#9b59b6')
    axes[2].set_xlabel('Days Elapsed')
    axes[2].set_ylabel('Clinical Episodes')
    axes[2].set_title('Absolute Clinical Episodes (Monthly)', fontweight='bold', fontsize=12)
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

# Main execution
# Set database path here (path to SQLite DB)
db_path = "monthly_data_100.db"  

if not db_path:
    raise ValueError("Please set variable 'db_path' to the SQLite DB path before running this cell.")

try:
    print(f"📊 Connecting to database: {db_path}")
    conn = connect_to_db(db_path)

    print("📈 Extracting data for last 12 months...")
    df = get_last_year_data(conn)

    if df.empty:
        print("❌ No data found in the database")
    else:
        print(f"✅ Found {len(df)} months of data")
        summary = get_summary_statistics(df)
        print_monthly_data(df)
        print_summary(summary)

        # Always create/show plot
        plot_data(df)

    conn.close()
    print("\n✅ Analysis complete!")

except FileNotFoundError as e:
    print(f"❌ Error: {e}")
except sqlite3.Error as e:
    print(f"❌ Database error: {e}")
except Exception as e:
    print(f"❌ Unexpected error: {e}")
    import traceback
    traceback.print_exc()
