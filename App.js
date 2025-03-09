import React, { useState } from "react";
import axios from "axios";
import "bootstrap/dist/css/bootstrap.min.css";
import "./App.css";

// ✅ Define Flask Backend URL
const API_BASE_URL = "https://chem-bot2.onrender.com";

function App() {
    const [smiles, setSmiles] = useState("");
    const [response, setResponse] = useState(null);
    const [imageUrl, setImageUrl] = useState(null);
    const [generatedMolecules, setGeneratedMolecules] = useState([]);
    const [functionalGroups, setFunctionalGroups] = useState([]);
    const [errorMessage, setErrorMessage] = useState("");
    const [loadingState, setLoadingState] = useState({
        analyzing: false,
        generating: false,
        detecting: false
    });

    // ✅ Handle Molecule Analysis & 2D Visualization (With IUPAC Naming)
    const handleSubmit = async () => {
        setErrorMessage("");
        setResponse(null);
        setImageUrl(null);

        if (!smiles.trim()) {
            setErrorMessage("❌ Please enter a valid SMILES string.");
            return;
        }

        try {
            setLoadingState((prev) => ({ ...prev, analyzing: true }));

            const res = await axios.post(`${API_BASE_URL}/molecule-info`, { smiles });
            setResponse(res.data);

            const imageRes = await axios.post(
                `${API_BASE_URL}/molecule-image`,
                { smiles },
                { responseType: "blob" }
            );
            setImageUrl(URL.createObjectURL(imageRes.data));
        } catch (error) {
            setErrorMessage("❌ Invalid molecule or server error.");
        } finally {
            setLoadingState((prev) => ({ ...prev, analyzing: false }));
        }
    };

    // ✅ Handle AI-Based Molecule Generation
    const handleGenerateMolecule = async () => {
        setErrorMessage("");
        setGeneratedMolecules([]);

        if (!smiles.trim()) {
            setErrorMessage("❌ Please enter a valid base SMILES string.");
            return;
        }

        try {
            setLoadingState((prev) => ({ ...prev, generating: true }));

            const res = await axios.post(`${API_BASE_URL}/generate-molecule`, { smiles });

            if (res.data && res.data.generated_smiles) {
                setGeneratedMolecules(res.data.generated_smiles);
            } else {
                setErrorMessage("❌ AI failed to generate valid molecules.");
            }
        } catch (error) {
            setErrorMessage("❌ AI failed to generate molecules.");
        } finally {
            setLoadingState((prev) => ({ ...prev, generating: false }));
        }
    };

    // ✅ Handle Functional Group Detection
    const handleFunctionalGroups = async () => {
        setErrorMessage("");
        setFunctionalGroups([]);

        if (!smiles.trim()) {
            setErrorMessage("❌ Please enter a valid SMILES string.");
            return;
        }

        try {
            setLoadingState((prev) => ({ ...prev, detecting: true }));

            const res = await axios.post(`${API_BASE_URL}/functional-groups`, { smiles });

            if (res.data && res.data.functional_groups) {
                setFunctionalGroups(res.data.functional_groups);
            } else {
                setErrorMessage("❌ No functional groups detected.");
            }
        } catch (error) {
            setErrorMessage("❌ Functional group detection failed.");
        } finally {
            setLoadingState((prev) => ({ ...prev, detecting: false }));
        }
    };

    return (
        <div className="container mt-5">
            <h1 className="text-center text-primary">🧪 Chemistry Research Chatbot</h1>

            {/* Molecule Input Field */}
            <div className="input-group my-4">
                <input
                    type="text"
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="Enter SMILES (e.g., CCO)"
                    className="form-control"
                />
                <button onClick={handleSubmit} className="btn btn-success">Analyze</button>
            </div>

            {/* Buttons Row */}
            <div className="d-flex justify-content-center mb-3">
                <button onClick={handleGenerateMolecule} className="btn btn-primary mx-2">Generate Molecule</button>
                <button onClick={handleFunctionalGroups} className="btn btn-dark mx-2">Detect Functional Groups</button>
            </div>

            {/* Show Loading Indicator */}
            {Object.values(loadingState).some(Boolean) && <p className="text-center text-info">⏳ Processing...</p>}

            {/* Show Error Message */}
            {errorMessage && <div className="alert alert-danger mt-3">{errorMessage}</div>}

            {/* Molecule Analysis Results (Including IUPAC Naming) */}
            {response && (
                <div className="card p-3 mt-3">
                    <h3 className="text-dark">🔬 Molecular Data</h3>
                    <p><strong>SMILES:</strong> {response.smiles}</p>
                    <p><strong>IUPAC Name:</strong> {response.iupac_name || "N/A"}</p>
                    <p><strong>Molecular Weight:</strong> {response.molecular_weight} g/mol</p>
                    <p><strong>Number of Atoms:</strong> {response.num_atoms}</p>
                    <p><strong>Number of Bonds:</strong> {response.num_bonds}</p>
                    <p><strong>LogP:</strong> {response.logP}</p>
                </div>
            )}

            {/* 2D Molecule Visualization */}
            {imageUrl && (
                <div className="mt-4 text-center">
                    <h4>🖼️ Molecule Structure</h4>
                    <img src={imageUrl} alt="Molecule Structure" className="img-fluid border border-dark rounded mt-3" />
                </div>
            )}

            {/* AI-Generated Molecules */}
            {generatedMolecules.length > 0 && (
                <div className="card p-3 mt-3">
                    <h3>🤖 AI-Generated Molecules</h3>
                    <ul className="list-group">
                        {generatedMolecules.map((mol, i) => (
                            <li key={i} className="list-group-item">{mol}</li>
                        ))}
                    </ul>
                </div>
            )}

            {/* Functional Groups */}
            {functionalGroups.length > 0 && (
                <div className="card p-3 mt-3">
                    <h3>🔬 Functional Groups Detected</h3>
                    <ul className="list-group">
                        {functionalGroups.map((group, i) => (
                            <li key={i} className="list-group-item">{group}</li>
                        ))}
                    </ul>
                </div>
            )}
        </div>
    );
}

export default App;


