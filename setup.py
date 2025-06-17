from setuptools import find_packages, setup

setup(
    name="natcules",
    version="0.1.0",
    author="Dr. Aayush Gupta",
    author_email="your.email@example.com",
    description="AI-driven Ayurvedic drug discovery using cheminformatics and ML",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/aaayushg/Natcules",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    install_requires=[
        "rdkit",                   # Cheminformatics
        "scikit-learn",            # Machine learning
        "pandas",                  # Data handling
        "numpy",                   # Numerical operations
        "matplotlib",              # Plotting (optional but common)
        "seaborn",                 # Visualization (optional)
        "tqdm",                    # Progress bars
        "joblib",                  # Model persistence
        "xgboost",                 # For advanced ML if used
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Adjust if different
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)

