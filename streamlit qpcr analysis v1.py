import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
from datetime import datetime
import re

# Page config
st.set_page_config(page_title="qPCR Analysis Suite", layout="wide", initial_sidebar_state="expanded")

# Initialize session state
if 'data' not in st.session_state:
    st.session_state.data = None
if 'processed_data' not in st.session_state:
    st.session_state.processed_data = None
if 'sample_mapping' not in st.session_state:
    st.session_state.sample_mapping = {}
if 'analysis_templates' not in st.session_state:
    st.session_state.analysis_templates = {}
if 'current_graph' not in st.session_state:
    st.session_state.current_graph = None
if 'excluded_wells' not in st.session_state:
    st.session_state.excluded_wells = set()

# Efficacy database
EFFICACY_DB = {
    '항노화': {'genes': ['MMP-1'], 'cell': 'HS68', 'controls': 1},
    '탄력': {'genes': ['COL1A1', 'ELN', 'FBN-1'], 'cell': 'HS68', 'controls': 1},
    '재생': {'genes': ['Migration'], 'cell': 'HS68', 'controls': 1},
    '보습': {'genes': ['HAS3', 'AQP3'], 'cell': 'HaCaT', 'controls': 1},
    '장벽': {'genes': ['FLG', 'CLDN', 'IVL'], 'cell': 'HaCaT', 'controls': 1},
    '진정': {'genes': ['IL-1β', 'IL-6', 'TNFα'], 'cell': 'HaCaT', 'controls': 2},
    '가려움': {'genes': ['TSLP'], 'cell': 'HaCaT', 'controls': 1},
    '표피증식': {'genes': ['PCNA', 'KI67'], 'cell': 'HaCaT', 'controls': 1},
    '멜라닌억제': {'genes': ['MITF', 'TYR'], 'cell': 'B16F10', 'controls': 1},
    '지질억제': {'genes': ['SREBPa', 'SREBPc'], 'cell': 'SZ95', 'controls': 1},
}

class QPCRParser:
    """Smart parser for multiple qPCR data formats"""
    
    @staticmethod
    def detect_format(df):
        """Auto-detect data format by searching through the entire file"""
        # Search for 'Well' or 'Well Position' in any cell
        for idx, row in df.iterrows():
            row_str = ' '.join(row.astype(str).values)
            if 'Well Position' in row_str:
                return 'format1', idx
            elif row.iloc[0] == 'Well' and 'Sample Name' in row_str and 'Target Name' in row_str:
                # Check if it's format1 or format2 by looking at other columns
                if 'Cт' in row_str or 'ΔCт' in row_str:
                    return 'format2', idx
                else:
                    return 'format1', idx
        return 'unknown', 0
    
    @staticmethod
    def parse_format1(df, data_start):
        """Parse QuantStudio format"""
        try:
            # Use the detected start row
            df = df.iloc[data_start:].reset_index(drop=True)
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)
            
            # Handle different column name variations
            well_col = 'Well Position' if 'Well Position' in df.columns else 'Well'
            ct_col = None
            for col in ['CT', 'Ct', 'Cт']:
                if col in df.columns:
                    ct_col = col
                    break
            
            if ct_col is None:
                st.error("Could not find CT/Ct column")
                return None
            
            # Extract relevant columns with error handling
            parsed = pd.DataFrame({
                'Well': df[well_col] if well_col in df.columns else df.iloc[:, 0],
                'Sample': df['Sample Name'] if 'Sample Name' in df.columns else df.iloc[:, 2],
                'Target': df['Target Name'] if 'Target Name' in df.columns else df.iloc[:, 3],
                'CT': pd.to_numeric(df[ct_col], errors='coerce'),
                'CT_Mean': pd.to_numeric(df.get('Ct Mean', df.get('CT Mean', df.get('Cт Mean', [np.nan]*len(df)))), errors='coerce'),
                'CT_SD': pd.to_numeric(df.get('Ct SD', df.get('CT SD', df.get('Cт SD', [np.nan]*len(df)))), errors='coerce'),
                'RQ': pd.to_numeric(df.get('RQ', [np.nan]*len(df)), errors='coerce'),
                'Delta_CT': pd.to_numeric(df.get('Delta Ct', df.get('ΔCт', [np.nan]*len(df))), errors='coerce'),
                'Delta_CT_Mean': pd.to_numeric(df.get('Delta Ct Mean', df.get('ΔCт Mean', [np.nan]*len(df))), errors='coerce'),
            })
            
            # Remove rows where CT is NaN (non-data rows)
            parsed = parsed.dropna(subset=['CT'])
            
            # Remove rows where Sample or Target is empty
            parsed = parsed[parsed['Sample'].notna() & parsed['Target'].notna()]
            
            return parsed
            
        except Exception as e:
            st.error(f"Error parsing format1: {str(e)}")
            return None
    
    @staticmethod
    def parse_format2(df, data_start):
        """Parse StepOnePlus format"""
        try:
            # Use the detected start row
            df = df.iloc[data_start:].reset_index(drop=True)
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)
            
            parsed = pd.DataFrame({
                'Well': df['Well'],
                'Sample': df['Sample Name'],
                'Target': df['Target Name'],
                'CT': pd.to_numeric(df['Cт'], errors='coerce'),
                'CT_Mean': pd.to_numeric(df.get('Cт Mean', [np.nan]*len(df)), errors='coerce'),
                'CT_SD': pd.to_numeric(df.get('Cт SD', [np.nan]*len(df)), errors='coerce'),
                'RQ': pd.to_numeric(df.get('RQ', [np.nan]*len(df)), errors='coerce'),
                'Delta_CT': pd.to_numeric(df.get('ΔCт', [np.nan]*len(df)), errors='coerce'),
                'Delta_CT_Mean': pd.to_numeric(df.get('ΔCт Mean', [np.nan]*len(df)), errors='coerce'),
            })
            
            # Remove rows where CT is NaN
            parsed = parsed.dropna(subset=['CT'])
            
            # Remove rows where Sample or Target is empty
            parsed = parsed[parsed['Sample'].notna() & parsed['Target'].notna()]
            
            return parsed
            
        except Exception as e:
            st.error(f"Error parsing format2: {str(e)}")
            return None
    
    @staticmethod
    def parse(file):
        """Main parsing function with robust error handling"""
        try:
            # Try different encodings
            df = None
            for encoding in ['utf-8', 'latin-1', 'cp1252']:
                try:
                    df = pd.read_csv(file, encoding=encoding, low_memory=False, skip_blank_lines=False)
                    break
                except UnicodeDecodeError:
                    continue
            
            if df is None:
                st.error("Could not read file with any encoding")
                return None
            
            # Detect format
            format_type, data_start = QPCRParser.detect_format(df)
            
            if format_type == 'format1':
                return QPCRParser.parse_format1(df, data_start)
            elif format_type == 'format2':
                return QPCRParser.parse_format2(df, data_start)
            else:
                st.error(f"Unknown format. Could not find data table in file. Please ensure your CSV contains 'Well', 'Sample Name', and 'Target Name' columns.")
                
                # Show preview to help debug
                with st.expander("📋 Show file preview (first 20 rows)"):
                    st.dataframe(df.head(20))
                
                return None
                
        except Exception as e:
            st.error(f"Error parsing file: {str(e)}")
            st.exception(e)
            return None

class StatisticalAnalyzer:
    """Statistical analysis engine"""
    
    @staticmethod
    def calculate_ddct(data, reference_sample, reference_gene, excluded_wells):
        """Calculate ΔΔCt and relative expression"""
        results = []
        
        # Filter out excluded wells
        data = data[~data['Well'].isin(excluded_wells)].copy()
        
        # Get reference values
        ref_data = data[(data['Sample'] == reference_sample) & (data['Target'] == reference_gene)]
        if len(ref_data) == 0:
            st.error(f"Reference sample '{reference_sample}' or gene '{reference_gene}' not found")
            return None
        
        ref_ct_mean = ref_data['CT'].mean()
        
        # Calculate for each sample-target combination
        for sample in data['Sample'].unique():
            for target in data['Target'].unique():
                subset = data[(data['Sample'] == sample) & (data['Target'] == target)]
                if len(subset) == 0:
                    continue
                
                # Get housekeeping gene CT for this sample
                hk_data = data[(data['Sample'] == sample) & (data['Target'] == reference_gene)]
                if len(hk_data) == 0:
                    continue
                
                hk_ct_mean = hk_data['CT'].mean()
                target_ct_mean = subset['CT'].mean()
                target_ct_sd = subset['CT'].std()
                
                # ΔCt = Ct(target) - Ct(housekeeping)
                delta_ct = target_ct_mean - hk_ct_mean
                
                # ΔΔCt = ΔCt(sample) - ΔCt(reference)
                ref_hk_ct = ref_data['CT'].mean()
                ref_target_data = data[(data['Sample'] == reference_sample) & (data['Target'] == target)]
                if len(ref_target_data) > 0:
                    ref_target_ct = ref_target_data['CT'].mean()
                    ref_delta_ct = ref_target_ct - ref_hk_ct
                else:
                    ref_delta_ct = 0
                
                ddct = delta_ct - ref_delta_ct
                
                # Relative expression = 2^(-ΔΔCt)
                rel_expr = 2 ** (-ddct)
                
                results.append({
                    'Sample': sample,
                    'Target': target,
                    'CT_Mean': target_ct_mean,
                    'CT_SD': target_ct_sd,
                    'n_replicates': len(subset),
                    'Delta_CT': delta_ct,
                    'Delta_Delta_CT': ddct,
                    'Relative_Expression': rel_expr,
                    'SEM': target_ct_sd / np.sqrt(len(subset)) if len(subset) > 1 else 0
                })
        
        return pd.DataFrame(results)
    
    @staticmethod
    def perform_statistics(processed_data, control_sample):
        """Perform t-tests comparing each sample to control"""
        results = processed_data.copy()
        results['p_value'] = np.nan
        results['significance'] = ''
        
        for target in results['Target'].unique():
            control_expr = results[(results['Sample'] == control_sample) & 
                                  (results['Target'] == target)]['Relative_Expression'].values
            
            for sample in results['Sample'].unique():
                if sample == control_sample:
                    continue
                
                sample_expr = results[(results['Sample'] == sample) & 
                                     (results['Target'] == target)]['Relative_Expression'].values
                
                if len(control_expr) > 0 and len(sample_expr) > 0:
                    # Perform t-test (assuming we have raw replicates, otherwise use means)
                    t_stat, p_val = stats.ttest_ind(control_expr, sample_expr)
                    
                    idx = (results['Sample'] == sample) & (results['Target'] == target)
                    results.loc[idx, 'p_value'] = p_val
                    
                    # Significance stars
                    if p_val < 0.001:
                        results.loc[idx, 'significance'] = '***'
                    elif p_val < 0.01:
                        results.loc[idx, 'significance'] = '**'
                    elif p_val < 0.05:
                        results.loc[idx, 'significance'] = '*'
        
        return results

class GraphGenerator:
    """Publication-ready graph generator"""
    
    @staticmethod
    def create_bar_chart(data, title, ylabel, show_error_bars=True):
        """Create grouped bar chart with error bars"""
        fig = go.Figure()
        
        samples = data['Sample'].unique()
        targets = data['Target'].unique()
        
        for target in targets:
            target_data = data[data['Target'] == target]
            
            fig.add_trace(go.Bar(
                name=target,
                x=target_data['Sample'],
                y=target_data['Relative_Expression'],
                error_y=dict(
                    type='data',
                    array=target_data['SEM'] * 1.96 if show_error_bars else None,  # 95% CI
                    visible=show_error_bars
                ),
                text=target_data['significance'],
                textposition='outside',
                textfont=dict(size=16)
            ))
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=20)),
            xaxis_title="Sample",
            yaxis_title=ylabel,
            barmode='group',
            template='plotly_white',
            font=dict(size=14),
            height=600,
            showlegend=True,
            legend=dict(font=dict(size=14))
        )
        
        return fig
    
    @staticmethod
    def create_heatmap(data):
        """Create heatmap for multi-gene analysis"""
        pivot = data.pivot(index='Target', columns='Sample', values='Relative_Expression')
        
        fig = go.Figure(data=go.Heatmap(
            z=pivot.values,
            x=pivot.columns,
            y=pivot.index,
            colorscale='RdYlGn',
            text=np.round(pivot.values, 2),
            texttemplate='%{text}',
            textfont={"size": 12},
            colorbar=dict(title="Fold Change")
        ))
        
        fig.update_layout(
            title="Gene Expression Heatmap",
            xaxis_title="Sample",
            yaxis_title="Gene",
            height=400 + len(pivot.index) * 30,
            font=dict(size=14)
        )
        
        return fig

def export_to_excel(raw_data, processed_data, analysis_params):
    """Export comprehensive Excel file with all calculations"""
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # Sheet 1: Analysis parameters
        params_df = pd.DataFrame([analysis_params])
        params_df.to_excel(writer, sheet_name='Analysis_Parameters', index=False)
        
        # Sheet 2: Raw data
        raw_data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        # Sheet 3: Processed data with calculations
        processed_data.to_excel(writer, sheet_name='Calculations', index=False)
        
        # Sheet 4: Summary statistics
        summary = processed_data.groupby('Target').agg({
            'Relative_Expression': ['mean', 'std', 'sem'],
            'p_value': 'min'
        }).round(4)
        summary.to_excel(writer, sheet_name='Summary')
    
    return output.getvalue()

# ==================== STREAMLIT UI ====================

st.title("🧬 qPCR Analysis Automation Suite")
st.markdown("**Automated analysis from raw data to publication-ready graphs**")

# Sidebar - Chat Assistant
with st.sidebar:
    st.header("💬 Quick Commands")
    st.markdown("""
    **Tips for faster analysis:**
    - Upload multiple files at once
    - Use templates for repeated assays
    - Flag outliers instead of deleting
    - Save templates for team reuse
    """)
    
    # Template management
    st.subheader("📋 Analysis Templates")
    template_name = st.text_input("Save current analysis as:")
    if st.button("💾 Save Template") and st.session_state.sample_mapping:
        st.session_state.analysis_templates[template_name] = {
            'mapping': st.session_state.sample_mapping.copy(),
            'timestamp': datetime.now().isoformat()
        }
        st.success(f"Template '{template_name}' saved!")
    
    if st.session_state.analysis_templates:
        load_template = st.selectbox("Load template:", [""] + list(st.session_state.analysis_templates.keys()))
        if load_template:
            st.session_state.sample_mapping = st.session_state.analysis_templates[load_template]['mapping']
            st.info(f"Template '{load_template}' loaded!")

# Main content tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📁 Upload", "🗺️ Sample Mapping", "🔬 Analysis", "📊 Visualization", "📤 Export"])

# TAB 1: File Upload
with tab1:
    st.header("Step 1: Upload qPCR Data")
    
    uploaded_files = st.file_uploader(
        "Upload CSV file(s) from QuantStudio or StepOnePlus",
        type=['csv'],
        accept_multiple_files=True
    )
    
    if uploaded_files:
        all_data = []
        for file in uploaded_files:
            parsed = QPCRParser.parse(file)
            if parsed is not None:
                parsed['Source_File'] = file.name
                all_data.append(parsed)
                st.success(f"✅ Parsed {file.name}: {len(parsed)} wells")
        
        if all_data:
            st.session_state.data = pd.concat(all_data, ignore_index=True)
            
            # Display preview
            st.subheader("Data Preview")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Wells", len(st.session_state.data))
            col2.metric("Samples", st.session_state.data['Sample'].nunique())
            col3.metric("Genes", st.session_state.data['Target'].nunique())
            
            # Quality metrics
            st.subheader("Quality Metrics")
            qc_data = st.session_state.data.groupby('Target').agg({
                'CT': ['mean', 'std', 'min', 'max']
            }).round(2)
            st.dataframe(qc_data)
            
            # Flag outliers
            st.subheader("🚩 Outlier Detection")
            for target in st.session_state.data['Target'].unique():
                target_data = st.session_state.data[st.session_state.data['Target'] == target]
                q1 = target_data['CT'].quantile(0.25)
                q3 = target_data['CT'].quantile(0.75)
                iqr = q3 - q1
                outliers = target_data[(target_data['CT'] < q1 - 1.5*iqr) | (target_data['CT'] > q3 + 1.5*iqr)]
                
                if len(outliers) > 0:
                    st.warning(f"**{target}**: {len(outliers)} potential outliers detected")
                    with st.expander(f"View {target} outliers"):
                        st.dataframe(outliers[['Well', 'Sample', 'CT']])
                        if st.checkbox(f"Exclude these wells from {target} analysis", key=f"exclude_{target}"):
                            st.session_state.excluded_wells.update(outliers['Well'].tolist())
            
            st.dataframe(st.session_state.data.head(20))

# TAB 2: Sample Mapping
with tab2:
    st.header("Step 2: Map Sample Names to Conditions")
    
    if st.session_state.data is not None:
        st.markdown("**Assign biological meaning to sample labels**")
        
        # Auto-detect efficacy type
        detected_genes = st.session_state.data['Target'].unique()
        suggested_efficacy = None
        for efficacy, info in EFFICACY_DB.items():
            if any(gene in detected_genes for gene in info['genes']):
                suggested_efficacy = efficacy
                break
        
        col1, col2 = st.columns(2)
        with col1:
            efficacy_type = st.selectbox(
                "Efficacy Test Type",
                options=list(EFFICACY_DB.keys()),
                index=list(EFFICACY_DB.keys()).index(suggested_efficacy) if suggested_efficacy else 0
            )
        
        with col2:
            expected_genes = EFFICACY_DB[efficacy_type]['genes']
            st.info(f"Expected genes: {', '.join(expected_genes)}")
        
        # Sample mapping interface
        st.subheader("Sample Mapping")
        samples = sorted(st.session_state.data['Sample'].unique(), key=lambda x: (not str(x).isdigit(), str(x)))
        
        mapping_data = []
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample,
                    'group': 'Treatment',
                    'concentration': ''
                }
            
            col1, col2, col3, col4 = st.columns([2, 3, 2, 2])
            with col1:
                st.text(sample)
            with col2:
                condition = st.text_input(
                    "Condition Name",
                    value=st.session_state.sample_mapping[sample]['condition'],
                    key=f"cond_{sample}"
                )
                st.session_state.sample_mapping[sample]['condition'] = condition
            with col3:
                group = st.selectbox(
                    "Group",
                    ['Control', 'Positive Control', 'Treatment', 'Induced Control'],
                    index=['Control', 'Positive Control', 'Treatment', 'Induced Control'].index(
                        st.session_state.sample_mapping[sample]['group']
                    ),
                    key=f"group_{sample}"
                )
                st.session_state.sample_mapping[sample]['group'] = group
            with col4:
                conc = st.text_input(
                    "Concentration",
                    value=st.session_state.sample_mapping[sample]['concentration'],
                    key=f"conc_{sample}"
                )
                st.session_state.sample_mapping[sample]['concentration'] = conc
        
        # Show mapping summary
        st.subheader("Mapping Summary")
        mapping_df = pd.DataFrame([
            {'Original': k, **v} for k, v in st.session_state.sample_mapping.items()
        ])
        st.dataframe(mapping_df)
    else:
        st.warning("⚠️ Please upload data first")

# TAB 3: Analysis
with tab3:
    st.header("Step 3: Statistical Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping:
        # Select reference sample and gene
        col1, col2 = st.columns(2)
        
        with col1:
            control_samples = [k for k, v in st.session_state.sample_mapping.items() 
                             if v['group'] == 'Control']
            reference_sample = st.selectbox("Reference Sample (Negative Control)", control_samples)
        
        with col2:
            housekeeping_genes = [g for g in st.session_state.data['Target'].unique() 
                                 if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            reference_gene = st.selectbox("Housekeeping Gene", housekeeping_genes)
        
        if st.button("🔬 Run Analysis", type="primary"):
            with st.spinner("Calculating ΔΔCt and statistics..."):
                # Calculate ΔΔCt
                processed = StatisticalAnalyzer.calculate_ddct(
                    st.session_state.data,
                    reference_sample,
                    reference_gene,
                    st.session_state.excluded_wells
                )
                
                if processed is not None:
                    # Perform statistics
                    processed = StatisticalAnalyzer.perform_statistics(processed, reference_sample)
                    
                    # Apply sample mapping
                    processed['Condition'] = processed['Sample'].map(
                        lambda x: st.session_state.sample_mapping.get(x, {}).get('condition', x)
                    )
                    processed['Group'] = processed['Sample'].map(
                        lambda x: st.session_state.sample_mapping.get(x, {}).get('group', 'Unknown')
                    )
                    
                    st.session_state.processed_data = processed
                    st.success("✅ Analysis complete!")
        
        # Display results
        if st.session_state.processed_data is not None:
            st.subheader("Analysis Results")
            
            # Summary statistics
            col1, col2, col3 = st.columns(3)
            sig_count = (st.session_state.processed_data['p_value'] < 0.05).sum()
            col1.metric("Significant Results", f"{sig_count} / {len(st.session_state.processed_data)}")
            col2.metric("Mean Fold Change", 
                       f"{st.session_state.processed_data['Relative_Expression'].mean():.2f}x")
            col3.metric("Max Fold Change", 
                       f"{st.session_state.processed_data['Relative_Expression'].max():.2f}x")
            
            # Data table with filtering
            st.subheader("Detailed Results")
            show_filter = st.checkbox("Show only significant results (p < 0.05)")
            
            display_data = st.session_state.processed_data.copy()
            if show_filter:
                display_data = display_data[display_data['p_value'] < 0.05]
            
            st.dataframe(
                display_data.style.background_gradient(subset=['Relative_Expression'], cmap='RdYlGn')
                .format({'Relative_Expression': '{:.3f}', 'p_value': '{:.4f}', 'CT_Mean': '{:.2f}'})
            )
    else:
        st.warning("⚠️ Complete previous steps first")

# TAB 4: Visualization
with tab4:
    st.header("Step 4: Generate Publication-Ready Graphs")
    
    if st.session_state.processed_data is not None:
        # Initialize graph settings in session state
        if 'graph_settings' not in st.session_state:
            st.session_state.graph_settings = {
                'title': 'Relative Gene Expression',
                'ylabel': 'Fold Change (Relative to Control)',
                'xlabel': 'Sample',
                'show_error': True,
                'font_size': 14,
                'title_size': 20,
                'bar_width': 0.8,
                'color_scheme': 'plotly_white',
                'show_legend': True,
                'legend_position': 'top right',
                'show_grid': True,
                'bar_colors': None,
                'figure_width': 1000,
                'figure_height': 600,
                'show_significance': True,
                'sig_font_size': 16
            }
        
        # Graph type selection
        st.subheader("🎨 Graph Customization")
        
        tab_type, tab_style, tab_colors, tab_labels = st.tabs([
            "📊 Graph Type", "🎨 Style", "🌈 Colors", "📝 Labels"
        ])
        
        # TAB: Graph Type
        with tab_type:
            col1, col2 = st.columns(2)
            with col1:
                graph_type = st.selectbox("Graph Type", ['Bar Chart', 'Heatmap', 'Both'])
            with col2:
                if graph_type in ['Bar Chart', 'Both']:
                    chart_orientation = st.radio("Orientation", ['Vertical', 'Horizontal'])
        
        # TAB: Style Settings
        with tab_style:
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.session_state.graph_settings['font_size'] = st.slider(
                    "Axis Font Size", 10, 24, st.session_state.graph_settings['font_size']
                )
                st.session_state.graph_settings['title_size'] = st.slider(
                    "Title Font Size", 14, 32, st.session_state.graph_settings['title_size']
                )
                st.session_state.graph_settings['sig_font_size'] = st.slider(
                    "Significance Star Size", 10, 24, st.session_state.graph_settings['sig_font_size']
                )
            
            with col2:
                st.session_state.graph_settings['bar_width'] = st.slider(
                    "Bar Width", 0.1, 1.0, st.session_state.graph_settings['bar_width']
                )
                st.session_state.graph_settings['figure_width'] = st.slider(
                    "Figure Width (px)", 600, 1600, st.session_state.graph_settings['figure_width']
                )
                st.session_state.graph_settings['figure_height'] = st.slider(
                    "Figure Height (px)", 400, 1200, st.session_state.graph_settings['figure_height']
                )
            
            with col3:
                st.session_state.graph_settings['color_scheme'] = st.selectbox(
                    "Color Theme",
                    ['plotly_white', 'plotly', 'plotly_dark', 'ggplot2', 'seaborn', 
                     'simple_white', 'presentation', 'none'],
                    index=['plotly_white', 'plotly', 'plotly_dark', 'ggplot2', 'seaborn', 
                           'simple_white', 'presentation', 'none'].index(
                        st.session_state.graph_settings['color_scheme']
                    )
                )
                st.session_state.graph_settings['show_legend'] = st.checkbox(
                    "Show Legend", st.session_state.graph_settings['show_legend']
                )
                st.session_state.graph_settings['show_grid'] = st.checkbox(
                    "Show Grid", st.session_state.graph_settings['show_grid']
                )
                st.session_state.graph_settings['show_error'] = st.checkbox(
                    "Show Error Bars", st.session_state.graph_settings['show_error']
                )
                st.session_state.graph_settings['show_significance'] = st.checkbox(
                    "Show Significance Stars", st.session_state.graph_settings['show_significance']
                )
        
        # TAB: Color Customization
        with tab_colors:
            st.markdown("**Custom Colors for Each Gene/Target**")
            
            targets = st.session_state.processed_data['Target'].unique()
            if st.session_state.graph_settings['bar_colors'] is None:
                # Initialize default colors
                default_colors = px.colors.qualitative.Plotly
                st.session_state.graph_settings['bar_colors'] = {
                    target: default_colors[i % len(default_colors)] 
                    for i, target in enumerate(targets)
                }
            
            cols = st.columns(min(3, len(targets)))
            for i, target in enumerate(targets):
                with cols[i % 3]:
                    st.session_state.graph_settings['bar_colors'][target] = st.color_picker(
                        f"{target}",
                        st.session_state.graph_settings['bar_colors'].get(target, '#636EFA')
                    )
            
            if st.button("🎨 Reset to Default Colors"):
                default_colors = px.colors.qualitative.Plotly
                st.session_state.graph_settings['bar_colors'] = {
                    target: default_colors[i % len(default_colors)] 
                    for i, target in enumerate(targets)
                }
                st.rerun()
        
        # TAB: Labels
        with tab_labels:
            col1, col2 = st.columns(2)
            
            with col1:
                st.session_state.graph_settings['title'] = st.text_input(
                    "Graph Title",
                    st.session_state.graph_settings['title']
                )
                st.session_state.graph_settings['xlabel'] = st.text_input(
                    "X-axis Label",
                    st.session_state.graph_settings['xlabel']
                )
                st.session_state.graph_settings['ylabel'] = st.text_input(
                    "Y-axis Label",
                    st.session_state.graph_settings['ylabel']
                )
            
            with col2:
                if st.session_state.graph_settings['show_legend']:
                    st.session_state.graph_settings['legend_position'] = st.selectbox(
                        "Legend Position",
                        ['top right', 'top left', 'bottom right', 'bottom left', 
                         'center right', 'center left'],
                        index=['top right', 'top left', 'bottom right', 'bottom left',
                               'center right', 'center left'].index(
                            st.session_state.graph_settings['legend_position']
                        )
                    )
        
        # Generate graph with all settings
        st.markdown("---")
        st.subheader("📊 Live Preview")
        
        if graph_type in ['Bar Chart', 'Both']:
            # Create bar chart with custom settings
            fig_bar = go.Figure()
            
            settings = st.session_state.graph_settings
            data = st.session_state.processed_data
            samples = data['Sample'].unique()
            targets = data['Target'].unique()
            
            for target in targets:
                target_data = data[data['Target'] == target]
                
                fig_bar.add_trace(go.Bar(
                    name=target,
                    x=target_data['Sample'] if chart_orientation == 'Vertical' else target_data['Relative_Expression'],
                    y=target_data['Relative_Expression'] if chart_orientation == 'Vertical' else target_data['Sample'],
                    orientation='v' if chart_orientation == 'Vertical' else 'h',
                    error_y=dict(
                        type='data',
                        array=target_data['SEM'] * 1.96 if settings['show_error'] and chart_orientation == 'Vertical' else None,
                        visible=settings['show_error'] and chart_orientation == 'Vertical'
                    ) if chart_orientation == 'Vertical' else None,
                    error_x=dict(
                        type='data',
                        array=target_data['SEM'] * 1.96 if settings['show_error'] and chart_orientation == 'Horizontal' else None,
                        visible=settings['show_error'] and chart_orientation == 'Horizontal'
                    ) if chart_orientation == 'Horizontal' else None,
                    text=target_data['significance'] if settings['show_significance'] else None,
                    textposition='outside',
                    textfont=dict(size=settings['sig_font_size']),
                    marker=dict(color=settings['bar_colors'].get(target, '#636EFA')),
                    width=settings['bar_width'] if chart_orientation == 'Vertical' else None
                ))
            
            # Legend position mapping
            legend_positions = {
                'top right': dict(x=1, y=1, xanchor='right', yanchor='top'),
                'top left': dict(x=0, y=1, xanchor='left', yanchor='top'),
                'bottom right': dict(x=1, y=0, xanchor='right', yanchor='bottom'),
                'bottom left': dict(x=0, y=0, xanchor='left', yanchor='bottom'),
                'center right': dict(x=1, y=0.5, xanchor='right', yanchor='middle'),
                'center left': dict(x=0, y=0.5, xanchor='left', yanchor='middle'),
            }
            
            fig_bar.update_layout(
                title=dict(text=settings['title'], font=dict(size=settings['title_size'])),
                xaxis_title=settings['xlabel'] if chart_orientation == 'Vertical' else settings['ylabel'],
                yaxis_title=settings['ylabel'] if chart_orientation == 'Vertical' else settings['xlabel'],
                barmode='group',
                template=settings['color_scheme'],
                font=dict(size=settings['font_size']),
                height=settings['figure_height'],
                width=settings['figure_width'],
                showlegend=settings['show_legend'],
                legend=legend_positions[settings['legend_position']] if settings['show_legend'] else dict(),
                xaxis=dict(showgrid=settings['show_grid']),
                yaxis=dict(showgrid=settings['show_grid'])
            )
            
            st.plotly_chart(fig_bar, use_container_width=True)
            st.session_state.current_graph = fig_bar
        
        if graph_type in ['Heatmap', 'Both']:
            st.subheader("Heatmap")
            
            # Create customizable heatmap
            pivot = st.session_state.processed_data.pivot(
                index='Target', columns='Sample', values='Relative_Expression'
            )
            
            fig_heat = go.Figure(data=go.Heatmap(
                z=pivot.values,
                x=pivot.columns,
                y=pivot.index,
                colorscale='RdYlGn',
                text=np.round(pivot.values, 2),
                texttemplate='%{text}',
                textfont={"size": st.session_state.graph_settings['font_size']},
                colorbar=dict(title="Fold Change", titlefont=dict(size=st.session_state.graph_settings['font_size']))
            ))
            
            fig_heat.update_layout(
                title=st.session_state.graph_settings['title'] + " (Heatmap)",
                xaxis_title=st.session_state.graph_settings['xlabel'],
                yaxis_title="Gene",
                height=400 + len(pivot.index) * 30,
                width=st.session_state.graph_settings['figure_width'],
                font=dict(size=st.session_state.graph_settings['font_size'])
            )
            
            st.plotly_chart(fig_heat, use_container_width=True)
        
        # Quick presets
        st.markdown("---")
        st.subheader("⚡ Quick Presets")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            if st.button("📄 Publication Style"):
                st.session_state.graph_settings.update({
                    'color_scheme': 'simple_white',
                    'font_size': 14,
                    'title_size': 18,
                    'show_grid': True,
                    'figure_width': 1000,
                    'figure_height': 600
                })
                st.rerun()
        
        with col2:
            if st.button("🎨 Presentation Style"):
                st.session_state.graph_settings.update({
                    'color_scheme': 'presentation',
                    'font_size': 18,
                    'title_size': 24,
                    'show_grid': False,
                    'figure_width': 1200,
                    'figure_height': 700
                })
                st.rerun()
        
        with col3:
            if st.button("🌙 Dark Mode"):
                st.session_state.graph_settings.update({
                    'color_scheme': 'plotly_dark',
                    'font_size': 14,
                    'title_size': 20,
                    'show_grid': True
                })
                st.rerun()
        
        with col4:
            if st.button("🔄 Reset All"):
                del st.session_state.graph_settings
                st.rerun()
    
    else:
        st.warning("⚠️ Run analysis first")

# TAB 5: Export
with tab5:
    st.header("Step 5: Export Results")
    
    if st.session_state.processed_data is not None:
        st.subheader("Download Options")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Excel export
            st.markdown("### 📊 Excel (Full Analysis)")
            analysis_params = {
                'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
                'Reference_Sample': reference_sample if 'reference_sample' in locals() else 'N/A',
                'Housekeeping_Gene': reference_gene if 'reference_gene' in locals() else 'N/A',
                'Excluded_Wells': len(st.session_state.excluded_wells)
            }
            
            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params
            )
            
            st.download_button(
                label="📥 Download Excel Report",
                data=excel_data,
                file_name=f"qPCR_analysis_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        with col2:
            # Graph export using Plotly's built-in HTML export (no Kaleido needed)
            st.markdown("### 📈 Graph Export")
            
            if st.session_state.current_graph:
                # Export format selection
                export_format = st.radio(
                    "Export Format",
                    ['Interactive HTML', 'Static Image (PNG)', 'Vector (SVG)'],
                    help="HTML is recommended - works in browsers and presentations"
                )
                
                if export_format == 'Interactive HTML':
                    # Export as interactive HTML (no Kaleido needed!)
                    html_buffer = io.StringIO()
                    st.session_state.current_graph.write_html(html_buffer)
                    html_data = html_buffer.getvalue()
                    
                    st.download_button(
                        label="📥 Download Interactive HTML",
                        data=html_data,
                        file_name=f"qPCR_graph_{datetime.now().strftime('%Y%m%d_%H%M')}.html",
                        mime="text/html",
                        help="Open in any browser - fully interactive and zoomable!"
                    )
                    
                    st.info("💡 HTML graphs work in PowerPoint, Google Slides, and all browsers!")
                
                elif export_format in ['Static Image (PNG)', 'Vector (SVG)']:
                    # Try to export as image (requires Kaleido)
                    try:
                        img_format = 'png' if export_format == 'Static Image (PNG)' else 'svg'
                        
                        img_bytes = st.session_state.current_graph.to_image(
                            format=img_format,
                            width=st.session_state.graph_settings['figure_width'],
                            height=st.session_state.graph_settings['figure_height'],
                            scale=2
                        )
                        
                        st.download_button(
                            label=f"📥 Download {img_format.upper()}",
                            data=img_bytes,
                            file_name=f"qPCR_graph_{datetime.now().strftime('%Y%m%d_%H%M')}.{img_format}",
                            mime=f"image/{img_format}"
                        )
                        
                    except Exception as e:
                        st.error("⚠️ Image export requires Kaleido (not available in cloud)")
                        st.info("**Workaround:** Download as HTML, then:")
                        st.markdown("""
                        1. **For PNG/JPG**: Open HTML in browser → Right-click graph → 'Save as image'
                        2. **For SVG**: Open HTML in browser → Click camera icon → Select 'svg' format
                        3. **Alternative**: Use screenshot tool for quick captures
                        """)
                        
                        # Fallback: Show JSON data for offline rendering
                        with st.expander("🔧 Advanced: Export as JSON (for offline rendering)"):
                            json_data = st.session_state.current_graph.to_json()
                            st.download_button(
                                label="📥 Download Graph JSON",
                                data=json_data,
                                file_name=f"qPCR_graph_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                                mime="application/json"
                            )
                            st.caption("Use this JSON with Plotly offline to render images locally")
        
        # CSV export for raw calculations
        st.markdown("---")
        st.subheader("📋 Additional Exports")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # Export processed data as CSV
            csv_buffer = io.StringIO()
            st.session_state.processed_data.to_csv(csv_buffer, index=False)
            
            st.download_button(
                label="📥 Download Calculations (CSV)",
                data=csv_buffer.getvalue(),
                file_name=f"qPCR_calculations_{datetime.now().strftime('%Y%m%d_%H%M')}.csv",
                mime="text/csv"
            )
        
        with col2:
            # JSON export for reproducibility
            repro_data = {
                'analysis_date': datetime.now().isoformat(),
                'sample_mapping': st.session_state.sample_mapping,
                'excluded_wells': list(st.session_state.excluded_wells),
                'graph_settings': st.session_state.graph_settings,
                'parameters': analysis_params if 'analysis_params' in locals() else {}
            }
            
            st.download_button(
                label="📥 Download Config (JSON)",
                data=json.dumps(repro_data, indent=2),
                file_name=f"qPCR_config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                mime="application/json"
            )
        
        with col3:
            # Export graph settings preset
            if st.button("💾 Save Graph Preset"):
                preset_name = st.text_input("Preset Name", "my_preset")
                if preset_name:
                    preset_data = json.dumps(st.session_state.graph_settings, indent=2)
                    st.download_button(
                        label=f"📥 Download '{preset_name}'",
                        data=preset_data,
                        file_name=f"graph_preset_{preset_name}.json",
                        mime="application/json"
                    )
        
        st.success("✅ All export options ready!")
        
        # Export tips
        with st.expander("💡 Export Tips & Best Practices"):
            st.markdown("""
            ### For Publications
            - Download **Excel** for peer review (shows all calculations)
            - Use **HTML** for interactive figures in supplementary materials
            - Convert HTML to PNG/SVG in browser for print journals
            
            ### For Presentations
            - **HTML** works directly in PowerPoint/Google Slides (drag & drop)
            - Interactive graphs impress audiences!
            - Or screenshot from browser for static slides
            
            ### For Patents/IP
            - Download **Excel** (complete audit trail)
            - Download **Config JSON** (reproducibility proof)
            - Save **HTML** (timestamped visual evidence)
            
            ### Browser Screenshot Method
            1. Download and open HTML file
            2. **Chrome/Edge**: Right-click graph → "Save image as..."
            3. **Firefox**: Right-click → "Take Screenshot" → Save
            4. **Safari**: Right-click → "Export as Image..."
            """)
    else:
        st.warning("⚠️ Complete analysis first")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>🧬 qPCR Analysis Suite v1.0 | Built for precision cosmetic science research</p>
    <p>Save templates • Flag outliers • Generate publication-ready graphs</p>
</div>
""", unsafe_allow_html=True)
