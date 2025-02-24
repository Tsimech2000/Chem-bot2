from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import random
from io import BytesIO

app = Flask(__name__)
CORS(app)

# ✅ Home Route
@app.route("/")
def home():
    return "Welcome to the Chemistry Research Chatbot!"

# ✅ Molecular Information Route
@app.route("/molecule-info", methods=["POST"])
def molecule_info():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "Invalid SMILES string"}), 400

        mol_Hs = Chem.AddHs(mol)

        properties = {
            "message": "Valid molecule",
            "smiles": smiles,
            "molecular_weight": round(Descriptors.MolWt(mol_Hs), 4),
            "num_atoms": mol_Hs.GetNumAtoms(),
            "num_bonds": mol_Hs.GetNumBonds(),
            "logP": round(Descriptors.MolLogP(mol_Hs), 4)
        }
        return jsonify(properties)

    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500

# ✅ Molecule 2D Image Generation Route
@app.route("/molecule-image", methods=["POST"])
def molecule_image():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "Invalid SMILES string"}), 400

        mol_Hs = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol_Hs)

        img_io = BytesIO()
        img = Draw.MolToImage(mol_Hs, size=(400, 400))
        img.save(img_io, format="PNG")
        img_io.seek(0)

        return send_file(img_io, mimetype="image/png")

    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500

# ✅ AI-Based Molecule Generation Route
@app.route("/generate-molecule", methods=["POST"])
def generate_molecule():
    """Generates new molecules by making small mutations to the input molecule."""
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No base SMILES input provided"}), 400

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "Invalid SMILES string"}), 400

        new_smiles = []
        for _ in range(5):  # Generate 5 new variations
            mutated_mol = Chem.RWMol(mol)

            atom_count = mutated_mol.GetNumAtoms()
            if atom_count > 1:
                random_atom = random.randint(0, atom_count - 1)
                mutated_mol.RemoveAtom(random_atom)
            AllChem.Compute2DCoords(mutated_mol)
            new_smiles.append(Chem.MolToSmiles(mutated_mol))

        valid_smiles = [smi for smi in new_smiles if Chem.MolFromSmiles(smi) is not None]

        if not valid_smiles:
            return jsonify({"error": "Failed to generate valid molecules"}), 400

        return jsonify({"generated_smiles": valid_smiles})

    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500

if __name__ == "__main__":
    app.run(debug=True)
























