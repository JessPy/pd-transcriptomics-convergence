#!/usr/bin/env python3
import os
import sys
import glob
import pandas as pd
import py4cytoscape as p4c
import time
import json

def create_and_export_network(de_file, out_root, gene_name, comparison_name):
    print(f"\n>>> Processing {gene_name} | Comparison: {comparison_name}")
    
    try:
        df = pd.read_csv(de_file)
    except Exception as e:
        print(f"  [ERROR] Could not read {de_file}: {e}")
        return

    # Filter significant genes (padj < 0.05)
    sig_df = df[df['padj'] < 0.05].copy()
    if sig_df.empty:
        print(f"  [SKIP] No significant genes found.")
        return

    # Prepare sets
    tasks = [
        ('Up', sig_df[sig_df['log2FoldChange'] > 0].sort_values('padj').head(100)),
        ('Down', sig_df[sig_df['log2FoldChange'] < 0].sort_values('padj').head(100)),
        ('Combined', pd.concat([
            sig_df[sig_df['log2FoldChange'] > 0].sort_values('padj').head(100),
            sig_df[sig_df['log2FoldChange'] < 0].sort_values('padj').head(100)
        ]))
    ]

    for direction, selected_genes in tasks:
        if selected_genes.empty:
            print(f"  [SKIP] No {direction} genes found.")
            continue
            
        target_dir = os.path.join(out_root, gene_name, comparison_name, direction)
        os.makedirs(target_dir, exist_ok=True)
        img_path = os.path.join(target_dir, f"network_{direction}_styled.png")

        gene_list = ",".join(selected_genes['symbol'].dropna().astype(str).tolist())
        network_title = f"{gene_name}_{comparison_name}_{direction}"

        try:
            print(f"  Querying STRING for {len(selected_genes)} {direction} genes...")
            
            # ── Execute STRING Query ──
            string_cmd = (
                f'string protein query query="{gene_list}" '
                f'species="Homo sapiens" '
                f'limit=20 '
                f'cutoff=0.8 '
                f'networkType="full STRING network"'
            )
            p4c.commands.commands_run(string_cmd)
            time.sleep(1) 
            
            # Ensure unique network name
            existing_networks = p4c.get_network_list()
            final_title = network_title
            counter = 1
            while final_title in existing_networks:
                final_title = f"{network_title}_{counter}"
                counter += 1
            p4c.rename_network(final_title)
            
            # ── Apply Visual Styling ──
            selected_genes['abs_log2FC'] = selected_genes['log2FoldChange'].abs()
            
            # Upload metadata to Nodes
            success_load = False
            for target_col in ['display name', 'query term', 'shared name', 'name']:
                try:
                    p4c.load_table_data(selected_genes[:], 
                                        data_key_column='symbol', 
                                        table_key_column=target_col)
                    print(f"    [INFO] Node data loaded successfully using: {target_col}")
                    success_load = True
                    break
                except:
                    continue
            
            # Define the style
            style_name = f"Style_{direction}_{gene_name}_{int(time.time())}"
            p4c.create_visual_style(style_name)

            def safe_map(func, *args, **kwargs):
                try:
                    func(*args, **kwargs)
                except Exception as ex:
                    print(f"    [WARNING] {func.__name__} failed: {ex}")

            print("    Applying mappings...")
            
            # 1. Node Color (Viridis)
            v_values = [-7, -5, -2, 0, 2, 5, 7]
            v_colors = [ "#503FCF", "#6170C4", "#85A2D2", "#C2DAD9", "#E28383", "#B94A4A",  "#AE2020"]
            safe_map(p4c.set_node_color_default, '#D3D3D3', style_name=style_name)
            safe_map(p4c.set_node_color_mapping, 'log2FoldChange', v_values, v_colors, mapping_type='c', style_name=style_name)
            
            # 2. Node Size
            safe_map(p4c.set_node_size_mapping, 'abs_log2FC', [1.0, 7.0], [50, 100], mapping_type='c', style_name=style_name)
            
            # 3. Node Labels
            safe_map(p4c.set_node_label_mapping, 'display name', style_name=style_name)

            print("    Setting defaults...")
            safe_map(p4c.set_visual_property_default, {'NODE_LABEL_FONT_SIZE': 17}, style_name=style_name)
            safe_map(p4c.set_visual_property_default, {'NODE_LABEL_COLOR': '#333333'}, style_name=style_name)
            safe_map(p4c.set_visual_property_default, {'EDGE_STROKE_UNSELECTED_PAINT': '#CCCCCC'}, style_name=style_name)

            p4c.set_visual_style(style_name)

            # ── Apply Layout ──
            try:
                p4c.layout_network(layout_name='yfiles-circular')
            except:
                p4c.layout_network(layout_name='circular')

            p4c.save_session(os.path.join(target_dir, f"network_{direction}.cys"))

            # ── Export ──
            csv_path = os.path.join(target_dir, f"genes_{direction}.csv")
            selected_genes.to_csv(csv_path, index=False)
            p4c.export_image(img_path, type='PNG', overwrite_file=True)
            print(f"  [SUCCESS] {direction} Image and Table exported to: {target_dir}")

            #######################
            # Exporta tabela STRING

            node_table = p4c.tables.get_table_columns(table='node')

            node_table.to_csv(os.path.join(target_dir, f"string_node_{direction}.csv"), index=False)
            print("Tabela de nós salva com sucesso!")

            edge_table = p4c.tables.get_table_columns(table='edge')
            edge_table.to_csv(os.path.join(target_dir, f"string_edge_{direction}.csv"), index=False)
            print("Tabela de interações salva com sucesso!")

        except Exception as e:
            print(f"  [ERROR] Failed to create/style {direction} network: {e}")

def main():
    root_dir = "/home/jess/Documents/estudos/pd-transcriptomics-convergence"
    out_root = os.path.join(root_dir, "data/cytoscape")
    
    try:
        p4c.cytoscape_ping()
        print("Connected to Cytoscape.")
    except Exception:
        print("CRITICAL ERROR: Cytoscape is not running.")
        sys.exit(1)

    de_files = glob.glob(os.path.join(root_dir, "data/de/**/results/DE_results_*.csv"), recursive=True)
    de_files = [f for f in de_files if 'genotypes_detailed' in f or 'GAB1' not in f]
    
    for f in sorted(de_files):
        parts = f.split(os.sep)
        if 'results' in parts:
            idx = parts.index('results')
            gene = parts[idx-1]
        else:
            gene = parts[-3]
            
        comp = parts[-1].replace("DE_results_", "").replace(".csv", "")
        print(f"\n=== Processing genes: {gene} ===")
        create_and_export_network(f, out_root, gene, comp)

if __name__ == "__main__":
    main()
