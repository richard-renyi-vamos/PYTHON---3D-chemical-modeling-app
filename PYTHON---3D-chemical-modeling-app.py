import sys
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QLineEdit, QPushButton, QLabel
from rdkit import Chem
from rdkit.Chem import AllChem

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
        self.resize(1000, 800)

        # Main Layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Control Panel
        controls = QHBoxLayout()
        self.input_field = QLineEdit()
        self.input_field.setPlaceholderText("Enter SMILES string (e.g., C1=CC=CC=C1 for Benzene)")
        self.btn_render = QPushButton("Render Molecule")
        self.btn_render.clicked.connect(self.update_molecule)
        
        controls.addWidget(QLabel("SMILES:"))
        controls.addWidget(self.input_field)
        controls.addWidget(self.btn_render)
        layout.addLayout(controls)

        # 3D Render Window
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)
        self.plotter.set_background("black")

    def generate_molecule_data(self, smiles):
        """Converts SMILES to 3D coordinates using RDKit."""
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add Hydrogens
        AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # Generate 3D coordinates
        return mol

    def update_molecule(self):
        smiles = self.input_field.text()
        try:
            mol = self.generate_molecule_data(smiles)
            self.render_molecule(mol)
        except Exception as e:
            print(f"Error: {e}")

    def render_molecule(self, mol):
        self.plotter.clear()
        conf = mol.GetConformer()

        # Render Atoms
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            color = ELEMENT_COLORS.get(symbol, "magenta")
            
            # Draw Sphere
            sphere = pv.Sphere(radius=0.3, center=(pos.x, pos.y, pos.z))
            self.plotter.add_mesh(sphere, color=color, smooth_shading=True)

        # Render Bonds
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            p1 = conf.GetAtomPosition(idx1)
            p2 = conf.GetAtomPosition(idx2)
            
            # Draw Line/Cylinder for bonds
            line = pv.Line((p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z))
            self.plotter.add_mesh(line, color="white", line_width=5)

        self.plotter.reset_camera()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MoleculeApp()
    window.show()
    sys.exit(app.exec_())
