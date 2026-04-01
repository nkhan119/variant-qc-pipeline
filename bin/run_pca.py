#!/usr/bin/env python3
import argparse
import pandas as pd
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pq", nargs="+", required=True, help="List of cohort parquets")
    parser.add_argument("--out-html", required=True)
    args = parser.parse_args()

    # 1. Load and Merge Data
    dfs = []
    for p in args.input_pq:
        try:
            temp_df = pd.read_parquet(p)
            # Standardize columns to uppercase
            temp_df.columns = [str(c).upper() for c in temp_df.columns]
            
            # If COHORT column is missing, infer it from filename (e.g., "Cohort_A")
            if 'COHORT' not in temp_df.columns:
                cohort_name = os.path.basename(p).split('_')[0]
                temp_df['COHORT'] = cohort_name
                
            dfs.append(temp_df)
            print(f"Successfully loaded: {p} (Cohort: {temp_df['COHORT'].iloc[0]})")
        except Exception as e:
            print(f"Error reading {p}: {e}")

    if not dfs:
        print("No valid parquet data found. Exiting.")
        sys.exit(1)

    df = pd.concat(dfs, ignore_index=True)

    # 2. Define Features for PCA
    features = ["SNV_COUNT", "TITV", "MEAN_DP", "HET_COUNT"]
    # Verify these exist in the dataframe
    available_features = [f for f in features if f in df.columns]
    
    if not available_features:
        print(f"Error: None of the expected features {features} were found in the data.")
        print(f"Available columns are: {list(df.columns)}")
        sys.exit(1)

    # 3. Pivot and Clean
    # Ensure SAMPLE also exists
    sample_col = 'SAMPLE' if 'SAMPLE' in df.columns else df.columns[0]
    
    pivot_df = df.set_index([sample_col, "COHORT"])[available_features].dropna()

    if pivot_df.empty or len(pivot_df) < 2:
        print("Insufficient data after dropping NAs for PCA.")
        sys.exit(1)

    # 4. PCA Calculation
    x = StandardScaler().fit_transform(pivot_df)
    pca = PCA(n_components=2)
    components = pca.fit_transform(x)
    
    var_explained = pca.explained_variance_ratio_ * 100

    # 5. Generate Plotly HTML
    labels = {str(i): f"PC {i+1} ({var_explained[i]:.1f}%)" for i in range(2)}
    fig = px.scatter(
        components, x=0, y=1, 
        color=pivot_df.index.get_level_values("COHORT"),
        hover_name=pivot_df.index.get_level_values(sample_col),
        labels=labels,
        title="Technical PCA: Variant Quality Stratification",
        template="plotly_white",
        symbol=pivot_df.index.get_level_values("COHORT")
    )

    fig.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))
    fig.update_layout(legend_title_text='Cohort')

    fig.write_html(args.out_html)
    print(f"PCA report saved to: {args.out_html}")

if __name__ == "__main__":
    main()
