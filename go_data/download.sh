#!/bin/bash
# Download GO reference data.
#
# The ontology and the annotation corpus MUST come from the same GO release --
# information content is computed against the annotations, and a mismatched DAG
# silently changes the IC values. Both files will download to go_data/<release>/.
#
# Usage: bash go_data/download.sh [release]   (default: 2025-07-22)

set -euo pipefail

RELEASE="${1:-2025-07-22}"
BASE="https://release.geneontology.org/${RELEASE}"

# Resolve relative to this script, so it works from any working directory.
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/${RELEASE}"
mkdir -p "$DEST"

# OBO FILE (GO DAG)
echo "==> go.obo (${RELEASE})"
curl -fL --retry 3 -o "${DEST}/go.obo" "${BASE}/ontology/go.obo"

# ANNOTATION FILE (UniProt GOA, IEA evidence removed)
echo "==> filtered_goa_uniprot_all_noiea.gaf (${RELEASE})"
curl -fL --retry 3 -o "${DEST}/filtered_goa_uniprot_all_noiea.gaf.gz" \
    "${BASE}/annotations/filtered_goa_uniprot_all_noiea.gaf.gz"
gunzip -f "${DEST}/filtered_goa_uniprot_all_noiea.gaf.gz"

echo "Done. Data in ${DEST}"
