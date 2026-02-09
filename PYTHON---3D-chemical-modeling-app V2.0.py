import sys
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                             QWidget, QLineEdit, QPushButton, QLabel, QMessageBox)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Atomic colors (CPK coloring convention)
ELEMENT_COLORS = {
    "C": "grey",
    "O": "red",
    "N": "blue",
    "H": "white",
    "S": "yellow",
    "P": "orange",
}

class MoleculeApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("3D Chemical Modeler")
        self.resize(1200, 800)

        # Main Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Control Panel
        controls = QHBoxLayout()
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("Enter SMILES (e.g., C1=CC=CC=C1 for Benzene)")
        
        self.btn_render = QPushButton("Render Molecule")
        self.btn_render.clicked.connect(self.update_molecule)
        
        # FEATURE: Add a clear button to reset the view
        self.btn_clear = QPushButton("Clear")
        self.btn_clear.clicked.connect(self.clear_view)
        
        controls.addWidget(QLabel("SMILES:"))
        controls.addWidget(self.input_field)
        controls.addWidget(self.btn_render)
        controls.addWidget(self.btn_clear)
        layout.addLayout(controls)

        # FEATURE: Info Panel for molecular data (Weight, Formula)
        self.info_label = QLabel("Enter a SMILES string to see properties.")
        self.info_label.setStyleSheet("font-weight: bold; color: #4CAF50; margin: 10px;")
        layout.addWidget(self.info_label)

        # 3D Render Window
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)
        self.plotter.set_background("black")

    def generate_molecule_data(self, smiles):
        """Converts SMILES to 3D coordinates using RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        mol = Chem.AddHs(mol)  # Essential for health/biological molecules to show H-bonding sites
        AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Generate 3D coordinates
        return mol

    def update_molecule(self):
        smiles = self.input_field.text()
        try:
            mol = self.generate_molecule_data(smiles)
            self.update_info_panel(mol) # Update the text info
            self.render_molecule(mol)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not render molecule: {e}")

    # FEATURE: Logic to calculate and display molecular properties
    def update_info_panel(self, mol):
        """Calculates molecular weight and formula for the user."""
        mw = Descriptors.MolWt(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        self.info_label.setText(f"Formula: {formula}  |  Molecular Weight: {mw:.2f} g/mol")

    # FEATURE: Clear function
    def clear_view(self):
        """Resets the plotter and the input fields."""
        self.plotter.clear()
        self.input_field.clear()
        self.info_label.setText("Enter a SMILES string to see properties.")

    def render_molecule(self, mol):
        self.plotter.clear()
        conf = mol.GetConformer()

        # Render Atoms
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            color = ELEMENT_COLORS.get(symbol, "magenta")
            
            # Spheres represent the electron cloud (VDW radius approximation)
            radius = 0.35 if symbol != "H" else 0.2
            sphere = pv.Sphere(radius=radius, center=(pos.x, pos.y, pos.z))
            self.plotter.add_mesh(sphere, color=color, smooth_shading=True)

        # FEATURE: Enhanced Bond Rendering
        # This checks if a bond is single, double, or triple and adjusts thickness
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            p1 = conf.GetAtomPosition(idx1)
            p2 = conf.GetAtomPosition(idx2)
            
            # Determine line width based on bond order
            b_type = bond.GetBondType()
            if b_type == Chem.rdchem.BondType.SINGLE:
                width = 4
            elif b_type == Chem.rdchem.BondType.DOUBLE:
                width = 10
            elif b_type == Chem.rdchem.BondType.TRIPLE:
                width = 16
            else: # Aromatic or other
                width = 7
            
            line = pv.Line((p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z))
            self.plotter.add_mesh(line, color="silver", line_width=width)

        self.plotter.reset_camera()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MoleculeApp()
    window.show()
    sys.exit(app.exec_())
