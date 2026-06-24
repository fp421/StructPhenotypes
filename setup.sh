#!/bin/bash
set -e

echo "Setting up StructPhenotypes environment..."

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing dependencies..."
pip install -r requirements.txt

# Create data directories
echo "Creating data directories..."
mkdir -p data/alphamissense
mkdir -p data/gnomad
mkdir -p data/alphafold
mkdir -p data/clinvar
mkdir -p outputs

# Optional: Set up gcloud for BigQuery access
if command -v gcloud &> /dev/null; then
    echo ""
    echo "gcloud CLI found. To use BigQuery features:"
    echo "  gcloud auth application-default login"
else
    echo ""
    echo "gcloud CLI not found. For BigQuery features, install:"
    echo "  https://cloud.google.com/sdk/docs/install"
fi

echo ""
echo "✓ Setup complete!"
echo ""
echo "Available commands:"
echo "  python structpheno/get_alphamissense.py GENE_NAME  # Download AlphaMissense predictions"
echo "  python structpheno/get_gnomad.py GENE_NAME         # Download gnomAD variants"
echo "  python structpheno/get_alpha_fold.py GENE_NAME     # Download AlphaFold structure"
echo ""
echo "Example:"
echo "  python structpheno/get_gnomad.py SCN2A --report"
