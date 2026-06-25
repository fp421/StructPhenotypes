#!/bin/bash
# Setup script to download AlphaMissense predictions from Google Cloud Storage

set -e

echo "Setting up AlphaMissense predictions..."

# Check if Google Cloud SDK is installed
if ! command -v gsutil &> /dev/null; then
    echo "❌ Google Cloud SDK not found."
    echo ""
    echo "Installing Google Cloud SDK via Homebrew..."
    brew install google-cloud-sdk
fi

# Authenticate if needed
if ! gcloud auth application-default print-access-token &> /dev/null; then
    echo "Authenticating with Google Cloud..."
    gcloud auth application-default login
fi

# Create output directory
mkdir -p data/alphamissense

# Download the AlphaMissense predictions (amino acid substitutions version)
echo "📥 Downloading AlphaMissense amino acid substitutions (~1.1GB)..."
echo "   This may take a few minutes..."
gsutil -m cp gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz data/alphamissense/

echo "✅ Successfully downloaded AlphaMissense data!"
echo ""
echo "You can now run:"
echo "  python structpheno/get_alphamissense.py SCN2A"
echo "  python structpheno/get_alphamissense.py TP53"
echo "  python structpheno/get_alphamissense.py [YOUR_GENE_NAME]"
echo ""
echo "Each run will extract all missense variants for your protein of interest."
