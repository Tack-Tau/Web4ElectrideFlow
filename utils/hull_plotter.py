"""
Convex hull plotter with Plotly for web visualization

Converts matplotlib-based convex hull plots to interactive Plotly plots
"""

import json
import re
from math import gcd
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import matplotlib.cm as cm
from scipy.spatial import Delaunay

import plotly.graph_objects as go

from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram


def formula_to_html(formula: str) -> str:
    """Convert chemical formula to HTML format with subscripts for web display.
    
    Handles both simple subscripts (Ca3 -> Ca<sub>3</sub>) and 
    parentheses subscripts (Ca3(AlP2)2 -> Ca<sub>3</sub>(AlP<sub>2</sub>)<sub>2</sub>).
    Ignores subscript 1.
    """
    # First, handle closing parenthesis followed by number: )2 -> )<sub>2</sub>
    def replace_paren(m):
        count = m.group(1)
        if count == '1':
            return ')'
        else:
            return f')<sub>{count}</sub>'
    
    result = re.sub(r'\)(\d+)', replace_paren, formula)
    
    # Then handle element followed by number: Ca2 -> Ca<sub>2</sub>
    def replace_elem(m):
        element = m.group(1)
        count = m.group(2)
        if count == '1':
            return element
        else:
            return f'{element}<sub>{count}</sub>'
    
    result = re.sub(r'([A-Z][a-z]?)(\d+)', replace_elem, result)
    
    return result


# Valence electrons for excess electron calculation
VALENCE_ELECTRONS = {
    'H': 1, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
    'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
    'B': 3, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3,
    'Sc': 3, 'Y': 3,
    'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
    'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,
    'O': 2, 'S': 2, 'Se': 2, 'Te': 2, 'Po': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'At': 1,
}

ELECTRONEGATIVE_ELEMENTS = {'N', 'P', 'As', 'Sb', 'Bi', 'O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I', 'At'}


def calculate_nexcess_boundaries_binary(
    elem_A: str, 
    elem_C: str, 
    n_excess_max: float = 4.0, 
    max_atoms: int = 20
) -> Tuple[Optional[float], Optional[float]]:
    """Calculate composition boundaries for N_excess region in binary system"""
    
    if elem_A not in VALENCE_ELECTRONS or elem_C not in VALENCE_ELECTRONS:
        return None, None
    
    val_A = VALENCE_ELECTRONS[elem_A]
    val_C = VALENCE_ELECTRONS[elem_C]
    
    valid_x_C = []
    
    for l in range(1, max_atoms):
        for n in range(1, max_atoms):
            if l + n > max_atoms:
                continue
            
            g = gcd(l, n)
            l_p, n_p = l // g, n // g
            
            if (l_p, n_p) != (l, n):
                continue
            
            n_excess = val_A * l_p - val_C * n_p
            
            if 0 < n_excess <= n_excess_max:
                x_C = n_p / (l_p + n_p)
                valid_x_C.append(x_C)
    
    if not valid_x_C:
        return None, None
    
    return min(valid_x_C), max(valid_x_C)


def generate_binary_hull_plotly(
    chemsys: str,
    electride_data: List[Dict],
    mp_phases: List[Dict],
    e_hull_max: float = 0.05
) -> Dict[str, Any]:
    """
    Generate binary convex hull plot using Plotly (automatic mode)
    
    Args:
        chemsys: Chemical system (e.g., "Ca-P")
        electride_data: List of electride structure data
        mp_phases: List of MP reference phases
        e_hull_max: Maximum energy above hull to show
        
    Returns:
        Dictionary with Plotly figure JSON
    """
    elements = sorted(chemsys.split('-'))
    if len(elements) != 2:
        return {'error': 'Binary system required'}
    
    elem_A, elem_C = elements
    
    # Deduplicate electrides by composition (keep lowest energy per composition)
    comp_groups = {}
    for e_data in electride_data:
        comp = Composition(e_data['composition'])
        comp_key = comp.reduced_formula
        
        if comp_key not in comp_groups:
            comp_groups[comp_key] = e_data
        else:
            if e_data['energy_per_atom'] < comp_groups[comp_key]['energy_per_atom']:
                comp_groups[comp_key] = e_data
    
    deduplicated_electrides = list(comp_groups.values())
    
    # Create ComputedEntry objects
    entries = []
    entry_metadata = {}
    
    for e_data in deduplicated_electrides:
        comp = Composition(e_data['composition'])
        n_atoms = comp.num_atoms
        total_energy = e_data['energy_per_atom'] * n_atoms
        
        entry = ComputedEntry(
            composition=comp,
            energy=total_energy,
            entry_id=e_data['structure_id']
        )
        entries.append(entry)
        entry_metadata[entry.entry_id] = {
            'is_electride': True,
            'formula': e_data['composition'],
            'full_dft_e_hull': e_data.get('full_dft_e_hull', 0),
            'pearson_symbol': e_data.get('pearson_symbol', ''),
        }
    
    for phase in mp_phases:
        comp = Composition(phase['composition'])
        phase_chemsys = '-'.join(sorted([str(e) for e in comp.elements]))
        
        if phase_chemsys == chemsys or len(comp.elements) == 1:
            if any(str(e) in elements for e in comp.elements):
                total_energy = phase['energy']
                entry_id = phase.get('entry_id', phase.get('mp_id', 'unknown'))
                
                entry = ComputedEntry(
                    composition=comp,
                    energy=total_energy,
                    entry_id=entry_id
                )
                entries.append(entry)
                entry_metadata[entry_id] = {'is_electride': False}
    
    if not entries:
        return {'error': 'No entries found'}
    
    mp_entries = [e for e in entries if not entry_metadata[e.entry_id]['is_electride']]
    pd_mp_only = PhaseDiagram(mp_entries)
    pd_all = PhaseDiagram(entries)
    
    fig = go.Figure()
    
    # Plot hull lines
    mp_hull_data = []
    full_hull_data = []
    
    for entry in entries:
        comp = entry.composition
        x = comp.get_atomic_fraction(elem_C)
        formation_energy = pd_all.get_form_energy_per_atom(entry)
        
        if entry in pd_all.stable_entries:
            full_hull_data.append({'x': x, 'y': formation_energy})
        if entry in pd_mp_only.stable_entries:
            mp_hull_data.append({'x': x, 'y': formation_energy})
    
    mp_hull_data.sort(key=lambda d: d['x'])
    full_hull_data.sort(key=lambda d: d['x'])
    
    if len(mp_hull_data) >= 2:
        fig.add_trace(go.Scatter(
            x=[d['x'] for d in mp_hull_data],
            y=[d['y'] for d in mp_hull_data],
            mode='lines',
            line=dict(color='gray', width=2, dash='dash'),
            name='MP-only hull',
            hoverinfo='skip',
            showlegend=True
        ))
    
    if len(full_hull_data) >= 2:
        fig.add_trace(go.Scatter(
            x=[d['x'] for d in full_hull_data],
            y=[d['y'] for d in full_hull_data],
            mode='lines',
            line=dict(color='black', width=3),
            name='Convex hull',
            hoverinfo='skip',
            showlegend=True
        ))
    
    # Categorize phases for plotting (matching local script styling)
    mp_stable = []
    new_stable = []
    new_meta = []
    
    for entry in entries:
        comp = entry.composition
        x = comp.get_atomic_fraction(elem_C)
        formation_energy = pd_all.get_form_energy_per_atom(entry)
        metadata = entry_metadata[entry.entry_id]
        is_electride = metadata['is_electride']
        
        # Use pre-calculated full_dft_e_hull from database for electrides
        if is_electride:
            e_above_hull = metadata.get('full_dft_e_hull', 0)
        else:
            # For MP phases, compute from phase diagram
            e_above_hull = pd_all.get_e_above_hull(entry)
        
        is_on_hull = e_above_hull < 0.001
        
        phase_data = {
            'x': x,
            'y': formation_energy,
            'e_hull': e_above_hull,
            'formula': comp.reduced_formula,
            'entry_id': entry.entry_id,
            'pearson_symbol': metadata.get('pearson_symbol', '') if is_electride else ''
        }
        
        if not is_electride and is_on_hull:
            mp_stable.append(phase_data)
        elif is_electride:
            if is_on_hull:
                new_stable.append(phase_data)
            elif e_above_hull <= e_hull_max:
                new_meta.append(phase_data)
    
    # Plot MP stable phases with tab20 colors (matching local script)
    tab20_colors = [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
        '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
        '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
        '#17becf', '#9edae5'
    ]
    
    for i, d in enumerate(mp_stable):
        color = tab20_colors[i % len(tab20_colors)]
        formula_html = formula_to_html(d['formula'])
        fig.add_trace(go.Scatter(
            x=[d['x']],
            y=[d['y']],
            mode='markers',
            marker=dict(size=15, color=color, symbol='circle', line=dict(color='black', width=1.5)),
            name=formula_html,
            hovertext=f"{d['formula']}<br>E_form = {d['y']:.4f} eV/atom",
            hoverinfo='text',
            showlegend=True
        ))
    
    # Plot stable electrides with different markers (matching local script)
    stable_markers = ['triangle-up', 'diamond', 'triangle-down', 'pentagon', 
                     'hexagon', 'x', 'hexagram', 'star', 'triangle-left', 'triangle-right']
    
    for i, d in enumerate(new_stable):
        marker = stable_markers[i % len(stable_markers)]
        pearson = d['pearson_symbol']
        formula_html = formula_to_html(d['formula'])
        
        # Create label: pearson-formula_html or just formula_html if no pearson
        if pearson:
            label = f'{pearson}-{formula_html}'
        else:
            label = formula_html
        
        fig.add_trace(go.Scatter(
            x=[d['x']],
            y=[d['y']],
            mode='markers+text',
            marker=dict(size=16, color='blue', symbol=marker, line=dict(color='red', width=2.5)),
            text=[label],
            textposition='top center',
            textfont=dict(size=10, color='blue', family='Arial Black'),
            name=label,
            hovertext=f"{d['entry_id']}<br>{d['formula']}<br>E_hull = {d['e_hull']:.4f} eV/atom",
            hoverinfo='text',
            showlegend=True
        ))
    
    # Plot metastable electrides (orange squares, matching local script)
    if new_meta:
        meta_x = [d['x'] for d in new_meta]
        meta_y = [d['y'] for d in new_meta]
        meta_text = [d['entry_id'] for d in new_meta]
        meta_hover = [f"{d['entry_id']}<br>{d['formula']}<br>E_hull = {d['e_hull']:.4f} eV/atom" 
                     for d in new_meta]
        
        fig.add_trace(go.Scatter(
            x=meta_x,
            y=meta_y,
            mode='markers',
            marker=dict(size=14, color='orange', symbol='square', line=dict(color='red', width=2)),
            name='Metastable electrides',
            hovertext=meta_hover,
            hoverinfo='text',
            showlegend=False
        ))
    
    # Configure layout (matching local script styling)
    title_text = f'{chemsys}'
    if new_meta:
        title_text = f'{chemsys} (E<sub>hull</sub> ≤ {e_hull_max:.2f} eV/atom)'
    
    fig.update_layout(
        title=dict(text=title_text, font=dict(size=20, family='Arial Black', color='black')),
        xaxis=dict(
            title=dict(text=f'Composition: x in {elem_A}<sub>1-x</sub>{elem_C}<sub>x</sub>', 
                      font=dict(size=16, family='Arial Black', color='black')),
            showgrid=True,
            gridcolor='rgba(200,200,200,0.3)',
            gridwidth=1,
            griddash='dash',
            range=[-0.05, 1.05],
            tickfont=dict(color='black', size=14),
            linecolor='black'
        ),
        yaxis=dict(
            title=dict(text='Formation Energy (eV/atom)', 
                      font=dict(size=16, family='Arial Black', color='black')),
            showgrid=True,
            gridcolor='rgba(200,200,200,0.3)',
            gridwidth=1,
            griddash='dash',
            tickfont=dict(color='black', size=14),
            linecolor='black'
        ),
        hovermode='closest',
        showlegend=True,
        legend=dict(
            x=0.98,
            y=0.02,
            xanchor='right',
            yanchor='bottom',
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='black',
            borderwidth=1,
            font=dict(size=12, color='black')
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family='Arial', size=14, color='black')
    )
    
    return {
        'data': [trace.to_plotly_json() for trace in fig.data],
        'layout': fig.layout.to_plotly_json()
    }


def generate_ternary_hull_plotly(
    chemsys: str,
    electride_data: List[Dict],
    mp_phases: List[Dict],
    e_hull_max: float = 0.05
) -> Dict[str, Any]:
    """
    Generate ternary convex hull plot using Plotly (automatic mode)
    
    Args:
        chemsys: Chemical system (e.g., "K-B-O")
        electride_data: List of electride structure data
        mp_phases: List of MP reference phases
        e_hull_max: Maximum energy above hull to show
        
    Returns:
        Dictionary with Plotly figure JSON
    """
    elements = sorted(chemsys.split('-'))
    if len(elements) != 3:
        return {'error': 'Ternary system required'}
    
    # Deduplicate electrides by composition (keep lowest energy)
    comp_groups = {}
    for data in electride_data:
        comp = Composition(data['composition'])
        comp_key = comp.reduced_formula
        
        if comp_key not in comp_groups:
            comp_groups[comp_key] = data
        else:
            if data['energy_per_atom'] < comp_groups[comp_key]['energy_per_atom']:
                comp_groups[comp_key] = data
    
    deduplicated_electrides = list(comp_groups.values())
    
    entries = []
    entry_metadata = {}
    
    for data in deduplicated_electrides:
        comp = Composition(data['composition'])
        n_atoms = comp.num_atoms
        total_energy = data['energy_per_atom'] * n_atoms
        
        entry = ComputedEntry(
            composition=comp,
            energy=total_energy,
            entry_id=data.get('structure_id', data['composition'])
        )
        entries.append(entry)
        entry_metadata[entry.entry_id] = {
            'is_electride': True,
            'formula': data['composition'],
            'structure_id': data.get('structure_id'),
            'full_dft_e_hull': data.get('full_dft_e_hull', 0),
            'pearson_symbol': data.get('pearson_symbol', ''),
        }
    
    for phase in mp_phases:
        try:
            phase_chemsys = phase.get('chemsys', '')
            comp_dict = phase['composition']
            comp = Composition(comp_dict)
            
            phase_elements = list(comp_dict.keys())
            comp_elements = set(phase_elements)
            system_elements = set(elements)
            
            if comp_elements.issubset(system_elements):
                entry = ComputedEntry(
                    composition=comp,
                    energy=phase['energy'],
                    entry_id=phase.get('entry_id', phase.get('mp_id', 'unknown'))
                )
                entries.append(entry)
                entry_metadata[entry.entry_id] = {
                    'is_electride': False,
                    'mp_id': phase.get('mp_id'),
                }
        except Exception as e:
            print(f"Warning: Skipping MP phase {phase.get('mp_id','unknown')}: {e}")
    
    if len(entries) == 0:
        return {'error': 'No entries found for system'}
    
    # Create phase diagrams
    pd_all = PhaseDiagram(entries)
    mp_only_entries = [e for e in entries if not entry_metadata[e.entry_id]['is_electride']]
    pd_mp_only = PhaseDiagram(mp_only_entries) if mp_only_entries else pd_all
    
    # Create ternary plot using Plotly
    fig = go.Figure()
    
    # Categorize phases (matching local script)
    mp_stable = []
    new_stable = []
    new_meta = []
    
    for entry in entries:
        try:
            comp = entry.composition
            fracs = [comp.get_atomic_fraction(el) for el in elements]
            a, b, c = fracs
            
            formation_energy = pd_all.get_form_energy_per_atom(entry)
            metadata = entry_metadata[entry.entry_id]
            is_electride = metadata['is_electride']
            
            # Use pre-calculated full_dft_e_hull from database for electrides
            if is_electride:
                e_above_hull = metadata.get('full_dft_e_hull', 0)
            else:
                # For MP phases, compute from phase diagram
                e_above_hull = pd_all.get_e_above_hull(entry)
            
            is_on_hull = e_above_hull < 0.001
            
            phase_data = {
                'a': a, 'b': b, 'c': c,
                'e_hull': e_above_hull,
                'e_form': formation_energy,
                'formula': comp.reduced_formula,
                'entry_id': entry.entry_id,
                'pearson_symbol': metadata.get('pearson_symbol', '') if is_electride else ''
            }
            
            if not is_electride and is_on_hull:
                mp_stable.append(phase_data)
            elif is_electride:
                if is_on_hull:
                    new_stable.append(phase_data)
                elif e_above_hull <= e_hull_max:
                    new_meta.append(phase_data)
        except Exception as e:
            print(f"Warning: Error processing entry {entry.entry_id}: {e}")
    
    # Draw hull boundaries: convex hull first (zorder=2), then MP-only hull (zorder=3) on top
    # Get stable purple color from rainbow colormap at E=0
    cmap = cm.rainbow
    stable_purple_rgba = cmap(0.0)
    stable_purple = '#{:02x}{:02x}{:02x}'.format(int(stable_purple_rgba[0]*255), 
                                                   int(stable_purple_rgba[1]*255), 
                                                   int(stable_purple_rgba[2]*255))
    
    # Draw full convex hull facets first (purple solid, zorder=2)
    stable_entries_on_hull = []
    for e in entries:
        try:
            if pd_all.get_e_above_hull(e) < 0.001:
                stable_entries_on_hull.append(e)
        except:
            pass  # Skip entries that can't be decomposed
    
    if len(stable_entries_on_hull) >= 3:
        
        stable_coords = []
        for entry in stable_entries_on_hull:
            comp = entry.composition
            a_coord = comp.get_atomic_fraction(elements[0])
            b_coord = comp.get_atomic_fraction(elements[1])
            c_coord = comp.get_atomic_fraction(elements[2])
            stable_coords.append([a_coord, b_coord, c_coord])
        
        stable_coords = np.array(stable_coords)
        
        try:
            tri = Delaunay(stable_coords[:, :2])
            drawn_edges = set()
            first_line = True
            
            for simplex in tri.simplices:
                for i in range(3):
                    p1_idx = simplex[i]
                    p2_idx = simplex[(i+1) % 3]
                    edge = tuple(sorted([p1_idx, p2_idx]))
                    
                    if edge not in drawn_edges:
                        drawn_edges.add(edge)
                        p1 = stable_coords[p1_idx]
                        p2 = stable_coords[p2_idx]
                        
                        fig.add_trace(go.Scatterternary(
                            a=[p1[0], p2[0]],
                            b=[p1[1], p2[1]],
                            c=[p1[2], p2[2]],
                            mode='lines',
                            line=dict(color=stable_purple, width=2.5),
                            showlegend=first_line,
                            name='Convex hull' if first_line else None,
                            hoverinfo='skip'
                        ))
                        first_line = False
        except Exception as e:
            print(f"Warning: Could not draw hull facets: {e}")
    
    # Draw MP-only hull facets on top (black dashed, zorder=3)
    mp_stable_entries = []
    for e in entries:
        try:
            if pd_mp_only.get_e_above_hull(e) < 0.001:
                mp_stable_entries.append(e)
        except:
            pass  # Skip entries that can't be decomposed
    
    if len(mp_stable_entries) >= 3:
        
        mp_coords = []
        for entry in mp_stable_entries:
            comp = entry.composition
            a_coord = comp.get_atomic_fraction(elements[0])
            b_coord = comp.get_atomic_fraction(elements[1])
            c_coord = comp.get_atomic_fraction(elements[2])
            mp_coords.append([a_coord, b_coord, c_coord])
        
        mp_coords = np.array(mp_coords)
        
        try:
            tri = Delaunay(mp_coords[:, :2])
            drawn_edges = set()
            first_line = True
            
            for simplex in tri.simplices:
                for i in range(3):
                    p1_idx = simplex[i]
                    p2_idx = simplex[(i+1) % 3]
                    edge = tuple(sorted([p1_idx, p2_idx]))
                    
                    if edge not in drawn_edges:
                        drawn_edges.add(edge)
                        p1 = mp_coords[p1_idx]
                        p2 = mp_coords[p2_idx]
                        
                        fig.add_trace(go.Scatterternary(
                            a=[p1[0], p2[0]],
                            b=[p1[1], p2[1]],
                            c=[p1[2], p2[2]],
                            mode='lines',
                            line=dict(color='black', width=2.5, dash='dash'),
                            showlegend=first_line,
                            name='MP-only hull' if first_line else None,
                            hoverinfo='skip'
                        ))
                        first_line = False
        except Exception as e:
            print(f"Warning: Could not draw MP-only hull facets: {e}")
    
    # Plot MP stable phases with tab20 colors (matching local script)
    tab20_colors = [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
        '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
        '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
        '#17becf', '#9edae5'
    ]
    
    for i, d in enumerate(mp_stable):
        color = tab20_colors[i % len(tab20_colors)]
        formula_html = formula_to_html(d['formula'])
        fig.add_trace(go.Scatterternary(
            a=[d['a']], b=[d['b']], c=[d['c']],
            mode='markers',
            marker=dict(size=18, color=color, symbol='circle', line=dict(color='black', width=1.5)),
            name=formula_html,
            hovertext=f"{d['formula']}<br>E_form = {d['e_form']:.4f} eV/atom",
            hoverinfo='text',
            showlegend=True
        ))
    
    # Stable electrides with purple color (from rainbow colormap at E=0) and different markers
    stable_markers = ['triangle-up', 'diamond', 'triangle-down', 'pentagon',
                     'hexagon', 'x', 'hexagram', 'star', 'triangle-left', 'triangle-right']
    
    for i, d in enumerate(new_stable):
        marker = stable_markers[i % len(stable_markers)]
        pearson = d['pearson_symbol']
        formula_html = formula_to_html(d['formula'])
        
        # Create label: pearson-formula_html or just formula_html if no pearson
        if pearson:
            label = f'{pearson}-{formula_html}'
        else:
            label = formula_html
        
        fig.add_trace(go.Scatterternary(
            a=[d['a']], b=[d['b']], c=[d['c']],
            mode='markers+text',
            marker=dict(size=20, color=stable_purple, symbol=marker, line=dict(color='red', width=2.5)),
            text=[label],
            textposition='top center',
            textfont=dict(size=11, color=stable_purple, family='Arial Black'),
            name=label,
            hovertext=f"{d['entry_id']}<br>{d['formula']}<br>E_hull = {d['e_hull']:.4f} eV/atom",
            hoverinfo='text',
            showlegend=True
        ))
    
    # Metastable electrides with rainbow colormap and colorbar
    if new_meta:
        e_hull_values = [d['e_hull'] for d in new_meta]
        vmin, vmax = 0.0, max(e_hull_values)
        
        # Create custom colorscale matching matplotlib's rainbow
        # Sample matplotlib rainbow at various points
        n_samples = 256
        matplotlib_rainbow = []
        for i in range(n_samples):
            val = i / (n_samples - 1)
            rgba = cmap(val)
            rgb_hex = '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
            matplotlib_rainbow.append([val, rgb_hex])
        
        # Collect all metastable data for a single trace with colorbar
        meta_a = [d['a'] for d in new_meta]
        meta_b = [d['b'] for d in new_meta]
        meta_c = [d['c'] for d in new_meta]
        meta_colors = [d['e_hull'] for d in new_meta]
        meta_hover = [f"{d['entry_id']}<br>{d['formula']}<br>E_hull = {d['e_hull']:.4f} eV/atom" 
                     for d in new_meta]
        
        fig.add_trace(go.Scatterternary(
            a=meta_a,
            b=meta_b,
            c=meta_c,
            mode='markers',
            marker=dict(
                size=16,
                color=meta_colors,
                colorscale=matplotlib_rainbow,
                cmin=vmin,
                cmax=vmax,
                symbol='square',
                line=dict(color='red', width=2),
                colorbar=dict(
                    title=dict(text='Energy above hull<br>(eV/atom)', font=dict(size=14, family='Arial Black', color='black')),
                    thickness=15,
                    len=0.80,
                    x=1.02,
                    y=0.5,
                    yanchor='middle',
                    tickfont=dict(size=12, color='black'),
                    outlinewidth=1,
                    outlinecolor='black'
                )
            ),
            hovertext=meta_hover,
            hoverinfo='text',
            showlegend=False
        ))
    
    # Update layout (matching local script styling)
    title_text = f'{chemsys}'
    if new_meta:
        title_text = f'{chemsys} (E<sub>hull</sub> ≤ {e_hull_max:.2f} eV/atom)'
    
    fig.update_layout(
        ternary=dict(
            sum=1,
            aaxis=dict(
                title=dict(text=elements[0], font=dict(size=18, family='Arial Black', color='black')), 
                min=0,
                showticklabels=False,
                linecolor='black',
                linewidth=2
            ),
            baxis=dict(
                title=dict(text=elements[1], font=dict(size=18, family='Arial Black', color='black')), 
                min=0,
                showticklabels=False,
                linecolor='black',
                linewidth=2
            ),
            caxis=dict(
                title=dict(text=elements[2], font=dict(size=18, family='Arial Black', color='black')), 
                min=0,
                showticklabels=False,
                linecolor='black',
                linewidth=2
            ),
            bgcolor='white'
        ),
        title=dict(text=title_text, font=dict(size=20, family='Arial Black', color='black')),
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='black',
            borderwidth=1,
            font=dict(size=12, color='black')
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family='Arial', size=14, color='black'),
        width=900,
        height=900,
        margin=dict(l=80, r=120, t=100, b=80)
    )
    
    # Convert to JSON-serializable format
    fig_dict = fig.to_dict()
    
    return {
        'data': json.loads(json.dumps(fig_dict['data'], default=str)),
        'layout': json.loads(json.dumps(fig_dict['layout'], default=str))
    }
