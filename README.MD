# 🌿 Natcules

**AI-Driven Discovery of Natural Product-Based Therapeutics Inspired by Global Traditional Medicine**

Natcules is a computational platform that combines **AI, cheminformatics, and computational chemistry** to discover novel, safer, and more effective therapeutic compounds rooted in ancient systems like **Ayurveda, Traditional Chinese Medicine (TCM), Korean medicine, and Indigenous plant knowledge**.

We fragment known bioactive compounds from curated traditional medicine datasets and use AI to generate and prioritize new molecular structures with optimized drug-like properties, safety, and biological relevance.

---

## 🚀 Key Features

- Integrates datasets from **IMPPAT, TCMID, NPASS, NANPDB**, and other traditional sources
- Uses **fragmentation (RECAP, BRICS, Murcko)** and **combinatorial enumeration**
- Predicts **drug-likeness, bioactivity, ADMET, and synthetic feasibility**
- Leverages **ML models** (e.g., ChemBERTa, Mol2Vec, GNNs) for property prediction
- Supports **partner applications** in biotech, pharma, nutraceuticals, cosmetics, and functional foods

---

## 🧪 Example Applications

- Generate optimized leads from Ashwagandha or Artemisia-derived compounds
- Predict and rank molecule candidates for anti-inflammatory targets
- Suggest functional food additives with Ayurvedic origins

---

## 📁 Project Structure
natcules/
├── data/ # Raw and processed traditional medicine data
├── notebooks/ # Jupyter notebooks for exploration and modeling
├── src/ # Core modules (fragmentation, modeling, generation)
├── models/ # Saved or pretrained models
├── scripts/ # CLI utilities for batch processing or predictions
├── tests/ # Unit tests
├── requirements.txt # Python dependencies
├── LICENSE # Project license (e.g., MIT or Apache 2.0)
└── README.md # You are here!

## ⚙️ Installation
git clone https://github.com/your-username/natcules.git
cd natcules
pip install -r requirements.txt

## 🧠 Requirements
Python 3.8+
RDKit
scikit-learn
PyTorch or TensorFlow (depending on your ML model backend)
pandas, numpy, matplotlib

## 🤝 Contributions
We welcome collaborators working on AI in natural products, ethnopharmacology, or traditional medicine digitization. Please open an issue or reach out directly.

## 📜 License
This project is licensed under the MIT License.

## ✨ About the Founder
Aayush Gupta, PhD in Computational Chemistry
Industry experience at Schrödinger, Exscientia, and DeepCure
Passionate about bridging ancient healing systems with cutting-edge AI to build safer, smarter therapeutics.

