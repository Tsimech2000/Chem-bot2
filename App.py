from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, DataStructs
import random
from io import BytesIO
import requests

app = Flask(__name__)
CORS(app)

# ✅ Functional Group SMARTS Patterns
FUNCTIONAL_GROUPS = {
    "Hydroxyl (-OH)": "[OX2H]", 
    "Amine (-NH2, -NR2)": "[NX3;H2,H1;!$(NC=O)]",  
    "Carbonyl (C=O)": "[CX3]=[OX1]",  
    "Carboxyl (-COOH)": "C(=O)[OX2H1]",  
    "Alkene (C=C)": "C=C",  
    "Alkyne (C≡C)": "C#C",  
    "Aromatic Ring": "a",  
    "Sulfonyl (-SO2-)": "[SX4](=O)(=O)",  
    "Ester (-COOR)": "[CX3](=O)[OX2H0]",  
    "Ether (-O-)": "[OD2]([#6])[#6]"  
}

# ✅ Sample Database for Similarity Search
MOLECULE_DATABASE = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ethanol": "CCO",
    "Methanol": "CO"
}

# ✅ Compute Tanimoto Similarity for Similarity Search
def compute_similarity(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 and mol2:
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    return 0.0

# ✅ Detect Functional Groups
def detect_functional_groups(mol):
    detected_groups = []
    for name, smarts in FUNCTIONAL_GROUPS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            detected_groups.append(name)
    return detected_groups

# ✅ Fetch IUPAC Name from PubChem API
def get_iupac_name(smiles):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["IUPACName"]
    except:
        return "IUPAC name not found"

# ✅ Home Route
@app.route("/")
def home():
    return "Welcome to the Chemistry Research Chatbot!"

# ✅ Molecule Analysis Route
@app.route("/molecule-info", methods=["POST"])
def molecule_info():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    iupac_name = get_iupac_name(smiles)
    properties = {
        "message": "Valid molecule",
        "smiles": smiles,
        "iupac_name": iupac_name,
        "molecular_weight": round(Descriptors.MolWt(mol), 4),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "logP": round(Descriptors.MolLogP(mol), 4)
    }

    return jsonify(properties)

# ✅ Molecule 2D Image Generation Route (Newly Added Feature)
@app.route("/molecule-image", methods=["POST"])
def molecule_image():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    AllChem.Compute2DCoords(mol)
    img_io = BytesIO()
    img = Draw.MolToImage(mol, size=(400, 400))
    img.save(img_io, format="PNG")
    img_io.seek(0)

    return send_file(img_io, mimetype="image/png")

# ✅ AI-Based Molecule Generation Route
@app.route("/generate-molecule", methods=["POST"])
def generate_molecule():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    generated_smiles = []
    for _ in range(5):
        try:
            mutated_mol = Chem.RWMol(mol)
            if mutated_mol.GetNumAtoms() > 1:
                random_atom = random.randint(0, mutated_mol.GetNumAtoms() - 1)
                mutated_mol.RemoveAtom(random_atom)
                Chem.SanitizeMol(mutated_mol)
                new_smiles = Chem.MolToSmiles(mutated_mol)
                if Chem.MolFromSmiles(new_smiles):
                    generated_smiles.append(new_smiles)
        except:
            continue  # Skip invalid mutations

    if not generated_smiles:
        return jsonify({"error": "Failed to generate valid molecules"}), 400

    return jsonify({"generated_smiles": generated_smiles})

# ✅ Functional Group Detection Route
@app.route("/functional-groups", methods=["POST"])
def functional_groups():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    groups = detect_functional_groups(mol)
    return jsonify({"smiles": smiles, "functional_groups": groups})

# ✅ Find Similar Molecules Route
@app.route("/similar-molecules", methods=["POST"])
def find_similar_molecules():
    data = request.json
    smiles = data.get("smiles", "").strip()

    if not smiles:
        return jsonify({"error": "No SMILES input provided"}), 400

    try:
        similarities = {
            name: compute_similarity(smiles, db_smiles)
            for name, db_smiles in MOLECULE_DATABASE.items()
        }

        sorted_similarities = sorted(similarities.items(), key=lambda x: x[1], reverse=True)

        return jsonify({
            "input_smiles": smiles,
            "similar_molecules": [
                {"name": name, "similarity": round(score, 3)}
                for name, score in sorted_similarities if score > 0
            ]
        })

    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500

if __name__ == "__main__":
    app.run(debug=True)
