#!/usr/bin/env python3
"""
Inorganic Electride Database Web Application

Main Flask application for hosting the electride database on Render.com.
Based on the electride_flow workflow: https://github.com/Tack-Tau/electride_flow
"""

import os
import sys
import json
import csv
import re
import traceback
from pathlib import Path
from io import StringIO

# Mock tkinter to avoid GUI dependencies on headless server
try:
    import tkinter
except ImportError:
    import types
    sys.modules['tkinter'] = types.ModuleType('tkinter')
    sys.modules['_tkinter'] = types.ModuleType('_tkinter')

# Add the bundled ASE to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'ase_root'))

from flask import render_template, jsonify, redirect, send_file
from jinja2 import FileSystemLoader

from ase.db.app import DBApp
from ase.db import connect
from ase.db.core import get_key_descriptions
from ase.db.project import DatabaseProject
from ase.io import write

from pymatgen.core import Composition

from utils.hull_plotter import (
    generate_binary_hull_plotly,
    generate_ternary_hull_plotly
)


class FilteredDatabase:
    """Wrapper to filter duplicate structure_ids from database queries and handle chemsys wildcards"""
    
    def __init__(self, db, valid_structure_ids):
        self.db = db
        self.valid_structure_ids = valid_structure_ids
        self._build_valid_rowids()
    
    def _build_valid_rowids(self):
        """Build set of valid row IDs - keep only FIRST occurrence of each structure_id"""
        self.valid_rowids = set()
        seen_structure_ids = set()
        
        for row in self.db.select():
            sid = row.get('structure_id', f'row_{row.id}')
            # Only keep first occurrence of each structure_id
            if sid in self.valid_structure_ids and sid not in seen_structure_ids:
                self.valid_rowids.add(row.id)
                seen_structure_ids.add(sid)
    
    def _match_chemsys_wildcard(self, chemsys_value, pattern):
        """
        Check if chemsys matches wildcard pattern.
        
        Supports patterns like:
        - '*' - match all
        - 'Ca-*' - binary with Ca as first element
        - '*-Ca' - binary with Ca as second element
        - 'B-*-O' - ternary with B first and O last
        - '*-K-O' - ternary with K second and O last
        - 'B-K-*' - ternary with B first and K second
        """
        if '*' not in pattern:
            return chemsys_value == pattern
        
        # Match all
        if pattern == '*':
            return True
        
        # Split both chemsys and pattern by '-'
        chemsys_parts = chemsys_value.split('-')
        pattern_parts = pattern.split('-')
        
        # Must have same number of elements
        if len(chemsys_parts) != len(pattern_parts):
            return False
        
        # Check each position
        for chem_part, pat_part in zip(chemsys_parts, pattern_parts):
            if pat_part != '*' and chem_part != pat_part:
                return False
        
        return True
    
    def select(self, *args, **kwargs):
        """Filter select to only return valid structure_ids and handle chemsys wildcards"""
        # Check if query contains chemsys wildcard
        query_str = args[0] if args else kwargs.get('query', '')
        chemsys_wildcard = None
        
        if isinstance(query_str, str) and 'chemsys=' in query_str:
            # Extract chemsys value from query
            match = re.search(r'chemsys=([^\s,;]+)', query_str)
            if match:
                chemsys_value = match.group(1)
                if '*' in chemsys_value:
                    chemsys_wildcard = chemsys_value
                    # Remove chemsys from query since we'll filter manually
                    query_str = re.sub(r'chemsys=[^\s,;]+[,;]?', '', query_str).strip()
                    if query_str.endswith(',') or query_str.endswith(';'):
                        query_str = query_str[:-1].strip()
                    # Update args or kwargs
                    if args:
                        args = (query_str,) + args[1:]
                    else:
                        kwargs['query'] = query_str
        
        # If we have a wildcard, we need to handle limit/offset manually
        if chemsys_wildcard:
            # Extract and remove limit/offset from kwargs
            original_limit = kwargs.pop('limit', None)
            original_offset = kwargs.pop('offset', 0)
            
            # Get ALL matching rows (no limit)
            matching_rows = []
            for row in self.db.select(*args, **kwargs):
                if row.id in self.valid_rowids:
                    row_chemsys = row.get('chemsys', '')
                    if self._match_chemsys_wildcard(row_chemsys, chemsys_wildcard):
                        matching_rows.append(row)
            
            # Apply offset and limit manually
            if original_offset:
                matching_rows = matching_rows[original_offset:]
            if original_limit:
                matching_rows = matching_rows[:original_limit]
            
            for row in matching_rows:
                yield row
        else:
            # No wildcard - normal filtering
            for row in self.db.select(*args, **kwargs):
                if row.id in self.valid_rowids:
                    yield row
    
    def count(self, *args, **kwargs):
        """Return count of valid structures"""
        if args or kwargs:
            # Count with filter
            return sum(1 for _ in self.select(*args, **kwargs))
        return len(self.valid_rowids)
    
    def get(self, *args, **kwargs):
        """Pass through get requests (single row access)"""
        return self.db.get(*args, **kwargs)
    
    def __getattr__(self, name):
        """Delegate other attributes to wrapped database"""
        return getattr(self.db, name)


class ElectrideDBApp:
    """Extended DBApp for electride database with custom routes"""
    
    def __init__(self):
        # Initialize base DBApp
        self.db_app = DBApp()
        self.flask = self.db_app.flask
        
        # CRITICAL: Configure template folder BEFORE any databases are added
        # This must point to the directory containing row.html, layout.html, etc.
        ase_db_path = Path(__file__).parent / 'ase_root' / 'ase' / 'db'
        template_path = ase_db_path / 'templates'
        self.flask.template_folder = str(template_path)
        
        # Also configure Jinja loader to look in the right place
        self.flask.jinja_loader = FileSystemLoader(str(template_path))
        
        # Configure static folder for ASE
        ase_static_path = ase_db_path / 'static'
        self.flask.static_folder = str(ase_static_path)
        self.flask.static_url_path = '/static'
        
        # Configure Flask app
        self.flask.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100MB max upload
        self.flask.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
        
        # Database connections
        self.databases = {}
        self.data_dir = Path(__file__).parent / 'data'
        
        # JSON data cache
        self.json_cache = {}
        
        # Structure ID whitelist (from CSV files to filter duplicates)
        self.valid_structure_ids = {}
        
        # Load databases
        self._load_databases()
        
        # Add custom routes
        self._add_custom_routes()
    
    def _load_databases(self):
        """Load all database files"""
        print("Loading databases...")
        
        # Load valid structure IDs from CSV files (to filter duplicates)
        self._load_csv_structure_ids()
        
        # Define custom columns to display (matching final_electrides.csv)
        # Note: Only show 'composition' (reduced formula), not 'formula' (full formula)
        custom_columns = [
            'id',
            'composition',
            'space_group_number',  # Maps to spacegroup in CSV
            'pearson_symbol',
            'dft_e_hull',
            'full_dft_e_hull',
            'mattersim_e_hull',
            'N_excess',
            'band_gap',
            'e0025',
            'e05',
            'e10',
            'band0',
            'band1',
            'electride',
        ]
        
        db_configs = [
            ('binary_final', 'binary/final_electrides.db'),
            ('ternary_final', 'ternary/final_electrides.db'),
        ]
        
        for name, rel_path in db_configs:
            db_path = self.data_dir / rel_path
            if db_path.exists():
                try:
                    raw_db = connect(str(db_path))
                    
                    # Wrap database with filter for valid structure_ids
                    system_type = 'binary' if 'binary' in name else 'ternary'
                    valid_ids = self.valid_structure_ids.get(system_type, set())
                    
                    if valid_ids:
                        db = FilteredDatabase(raw_db, valid_ids)
                        unique_count = db.count()
                    else:
                        db = raw_db
                        unique_count = raw_db.count()
                    
                    # Add project with custom columns
                    self._add_project_with_custom_columns(name, db, custom_columns)
                    self.databases[name] = db
                    print(f"  Loaded {name}: {raw_db.count()} rows ({unique_count} unique from CSV)")
                except Exception as e:
                    print(f"  ERROR loading {name}: {e}")
            else:
                print(f"  WARNING: {db_path} not found")
        
        # Load JSON data
        self._load_json_data()
    
    def _load_csv_structure_ids(self):
        """Load valid structure IDs from CSV files to filter duplicates"""
        csv_configs = [
            ('binary', 'binary/final_electrides.csv'),
            ('ternary', 'ternary/final_electrides.csv'),
        ]
        
        for system_type, rel_path in csv_configs:
            csv_path = self.data_dir / rel_path
            if csv_path.exists():
                try:
                    with open(csv_path, 'r') as f:
                        reader = csv.DictReader(f)
                        structure_ids = set()
                        for row in reader:
                            # CSV uses 'formula' column for structure_id
                            sid = row.get('formula', '').strip()
                            if sid:
                                structure_ids.add(sid)
                        self.valid_structure_ids[system_type] = structure_ids
                        print(f"  Loaded {len(structure_ids)} valid structure IDs from {rel_path}")
                except Exception as e:
                    print(f"  ERROR loading CSV {rel_path}: {e}")
                    self.valid_structure_ids[system_type] = set()
            else:
                print(f"  WARNING: CSV not found: {csv_path}")
                self.valid_structure_ids[system_type] = set()
    
    def _add_project_with_custom_columns(self, name, db, custom_columns):
        """Add database project with custom default columns"""
        
        # Create project with custom columns
        project = DatabaseProject(
            name=name,
            title=db.metadata.get('title', ''),
            key_descriptions=get_key_descriptions(),
            database=db,
            default_columns=custom_columns
        )
        
        # Add to Flask app
        self.db_app.projects[name] = project
    
    def _load_json_data(self):
        """Load JSON data files for each system"""
        print("Loading JSON data...")
        
        for system in ['binary', 'ternary']:
            system_dir = self.data_dir / system
            self.json_cache[system] = {}
            
            json_files = [
                'dft_stability_results.json',
                'mattersim_stability_results.json',
                'mp_vaspdft.json',
                'mp_mattersim.json',
            ]
            
            for json_file in json_files:
                json_path = system_dir / json_file
                if json_path.exists():
                    try:
                        with open(json_path, 'r') as f:
                            self.json_cache[system][json_file.replace('.json', '')] = json.load(f)
                        print(f"  Loaded {system}/{json_file}")
                    except Exception as e:
                        print(f"  ERROR loading {system}/{json_file}: {e}")
    
    def _add_custom_routes(self):
        """Add custom routes for electride-specific features"""
        
        @self.flask.route('/')
        def index():
            """Landing page with database statistics"""
            stats = self._get_statistics()
            return render_template('index.html', stats=stats)
        
        @self.flask.route('/browse/<project>')
        def browse(project):
            """Browse database with custom view"""
            if project not in self.databases:
                return f"Database '{project}' not found", 404
            return redirect(f'/{project}/')
        
        @self.flask.route('/structure/<project>/<int:structure_id>')
        def structure_detail(project, structure_id):
            """Detailed structure view with convex hull plot"""
            if project not in self.databases:
                return f"Database '{project}' not found", 404
            
            try:
                db = self.databases[project]
                row = db.get(id=structure_id)
                
                # Get system type and chemical system
                system_type = 'binary' if 'binary' in project else 'ternary'
                chemsys = self._get_chemical_system(row)
                
                # Get structure data
                structure_data = {
                    'id': structure_id,
                    'formula': row.get('structure_id', row.formula),
                    'composition': row.get('composition', row.formula),
                    'space_group': row.get('space_group_number', 'N/A'),
                    'dft_e_hull': row.get('dft_e_hull', None),
                    'mattersim_e_hull': row.get('mattersim_e_hull', None),
                    'full_dft_e_hull': row.get('full_dft_e_hull', None),
                    'band_gap': row.get('band_gap', None),
                    'N_excess': row.get('N_excess', None),
                    'e0025': row.get('e0025', None),
                    'e05': row.get('e05', None),
                    'e10': row.get('e10', None),
                    'band0': row.get('band0', None),
                    'band1': row.get('band1', None),
                    'electride': row.get('electride', None),
                }
                
                return render_template('structure.html',
                                      project=project,
                                      structure=structure_data,
                                      system_type=system_type,
                                      chemsys=chemsys)
            except Exception as e:
                return f"Error loading structure: {e}", 500
        
        @self.flask.route('/api/hull/<project>/<chemsys>')
        def get_hull_plot(project, chemsys):
            """Get convex hull plot data as JSON"""
            if project not in self.databases:
                return jsonify({'error': 'Database not found'}), 404
            
            try:
                system_type = 'binary' if 'binary' in project else 'ternary'
                plot_data = self._generate_hull_plot(project, system_type, chemsys)
                return jsonify(plot_data)
            except Exception as e:
                return jsonify({'error': str(e)}), 500
        
        @self.flask.route('/download/<project>/<format>')
        def download_database(project, format):
            """Download database in specified format"""
            if project not in self.databases:
                return "Database not found", 404
            
            if format not in ['db', 'csv', 'json']:
                return "Invalid format. Use 'db', 'csv', or 'json'", 400
            
            try:
                system_type = 'binary' if 'binary' in project else 'ternary'
                db_type = 'final' if 'final' in project else 'stable'
                
                if format == 'db':
                    db_path = self.data_dir / system_type / f'{db_type}_electrides.db'
                    return send_file(db_path, as_attachment=True,
                                   download_name=f'{project}.db')
                elif format == 'csv':
                    csv_path = self.data_dir / system_type / f'{db_type}_electrides.csv'
                    if csv_path.exists():
                        return send_file(csv_path, as_attachment=True,
                                       download_name=f'{project}.csv')
                    return "CSV file not found", 404
                elif format == 'json':
                    # Export database to JSON
                    db = self.databases[project]
                    rows = []
                    for row in db.select():
                        row_dict = dict(row.key_value_pairs)
                        row_dict['formula'] = row.formula
                        row_dict['id'] = row.id
                        rows.append(row_dict)
                    
                    return jsonify(rows)
            except Exception as e:
                return f"Error downloading: {e}", 500
        
        @self.flask.route('/atoms/<project>/<int:uid>/<format>')
        def get_atoms(project, uid, format):
            """Get structure file in specified format"""
            if project not in self.databases:
                return "Database not found", 404
            
            try:
                db = self.databases[project]
                row = db.get(id=uid)
                atoms = row.toatoms()
                
                if format == 'cif':
                    output = StringIO()
                    write(output, atoms, format='cif')
                    return output.getvalue(), 200, {'Content-Type': 'chemical/x-cif'}
                
                elif format == 'xyz':
                    output = StringIO()
                    write(output, atoms, format='xyz')
                    return output.getvalue(), 200, {'Content-Type': 'chemical/x-xyz'}
                
                elif format == 'json':
                    # Convert atoms to JSON-serializable dict
                    atoms_dict = {
                        'symbols': atoms.get_chemical_symbols(),
                        'positions': atoms.get_positions().tolist(),
                        'cell': atoms.get_cell().tolist(),
                        'pbc': atoms.get_pbc().tolist()
                    }
                    return jsonify(atoms_dict)
                
                else:
                    return "Invalid format. Use 'cif', 'xyz', or 'json'", 400
                    
            except Exception as e:
                return f"Error generating structure file: {e}", 500
        
        @self.flask.route('/about')
        def about():
            """About page with project information"""
            return render_template('about.html')
    
    def _get_statistics(self):
        """Calculate database statistics"""
        stats = {
            'total_structures': 0,
            'binary_final': 0,
            'binary_stable': 0,
            'ternary_final': 0,
            'ternary_stable': 0,
            'electrides': 0,
        }
        
        for name, db in self.databases.items():
            count = db.count()
            stats[name] = count
            stats['total_structures'] += count
            
            # Count confirmed electrides
            try:
                electride_count = db.count('electride=true')
                stats['electrides'] += electride_count
            except:
                pass
        
        return stats
    
    def _get_chemical_system(self, row):
        """Extract chemical system from database row"""
        try:
            atoms = row.toatoms()
            elements = sorted(set(atoms.get_chemical_symbols()))
            return '-'.join(elements)
        except:
            return 'unknown'
    
    def _generate_hull_plot(self, project, system_type, chemsys):
        """
        Generate convex hull plot for a SINGLE chemical system
        
        This method:
        1. Filters database to only include structures in the specified chemical system
        2. Generates ONE convex hull plot for that system
        3. Shows all competing phases with E_hull <= 0.05 eV/atom
        
        Args:
            project: Database project name (e.g., 'binary_final')
            system_type: 'binary' or 'ternary'
            chemsys: Chemical system string (e.g., 'Ca-P' or 'K-B-O')
        """
        try:
            db = self.databases[project]
            
            # Filter database to only include structures in this chemical system
            electride_data = []
            elements_set = set(chemsys.split('-'))
            
            for row in db.select():
                formula = row.formula
                comp = Composition(formula)
                row_elements = set(str(e) for e in comp.elements)
                
                if row_elements == elements_set:
                    electride_data.append({
                        'structure_id': row.get('structure_id', formula),
                        'composition': formula,
                        'energy_per_atom': row.get('vasp_energy_per_atom', row.get('mattersim_energy_per_atom', 0)),
                        'full_dft_e_hull': row.get('full_dft_e_hull', 0),
                        'pearson_symbol': row.get('pearson_symbol', ''),
                    })
            
            if not electride_data:
                return {
                    'project': project,
                    'system_type': system_type,
                    'chemsys': chemsys,
                    'data': [],
                    'layout': {},
                    'message': f'No structures found for chemical system {chemsys}'
                }
            
            json_key = 'binary' if 'binary' in project else 'ternary'
            mp_json_path = self.data_dir / json_key / 'mp_vaspdft.json'
            
            if not mp_json_path.exists():
                return {'error': f'MP reference phases file not found: {mp_json_path}'}
            
            with open(mp_json_path, 'r') as f:
                mp_phases = json.load(f)
            
            if system_type == 'binary':
                plot_json = generate_binary_hull_plotly(
                    chemsys=chemsys,
                    electride_data=electride_data,
                    mp_phases=mp_phases,
                    e_hull_max=0.05
                )
            else:
                plot_json = generate_ternary_hull_plotly(
                    chemsys=chemsys,
                    electride_data=electride_data,
                    mp_phases=mp_phases,
                    e_hull_max=0.05
                )
            
            return plot_json
            
        except Exception as e:
            return {
                'project': project,
                'system_type': system_type,
                'chemsys': chemsys,
                'error': str(e),
                'traceback': traceback.format_exc()
            }
    
    def run(self, host='0.0.0.0', port=None, debug=False):
        """Run the Flask application"""
        if port is None:
            port = int(os.environ.get('PORT', 5000))
        
        print(f"\nStarting Electride Database on http://{host}:{port}")
        print(f"Databases loaded: {len(self.databases)}")
        print(f"Total structures: {sum(db.count() for db in self.databases.values())}\n")
        
        self.flask.run(host=host, port=port, debug=debug, use_reloader=False)


def main():
    """Main entry point"""
    app = ElectrideDBApp()
    app.run()


if __name__ == '__main__':
    main()
