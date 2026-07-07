#!/usr/bin/env bash
# =============================================================================
# circRNA pipeline — Docker build, test, push, and Singularity deploy
#
# Usage:
#   On WSL2 (local):   bash containers/build_and_deploy.sh
#   On HPC server:     bash containers/build_and_deploy.sh --hpc
# =============================================================================

set -euo pipefail

DOCKERHUB_USER="choukaihsuan"
IMAGE_NAME="circrna-pipeline"
IMAGE_TAG="1.0.1"
FULL_IMAGE="${DOCKERHUB_USER}/${IMAGE_NAME}:${IMAGE_TAG}"
SIF_DIR="${HOME}/circRNA_agent/containers"
SIF_FILE="${SIF_DIR}/${IMAGE_NAME}_${IMAGE_TAG}.sif"

# ── HPC-only mode ─────────────────────────────────────────────────────────────
if [[ "${1:-}" == "--hpc" ]]; then
    echo "=== [HPC] Step 4: Pull Docker image and convert to Singularity .sif ==="
    mkdir -p "${SIF_DIR}"
    singularity pull "${SIF_FILE}" "docker://${FULL_IMAGE}"

    echo "=== [HPC] Step 5: Verify .sif ==="
    singularity exec "${SIF_FILE}" CIRIquant --version
    singularity exec "${SIF_FILE}" snakemake --version
    singularity exec "${SIF_FILE}" Rscript -e "library(edgeR); cat('edgeR OK\n')"

    echo ""
    echo "=== [HPC] Step 6: Run pipeline with Singularity ==="
    echo "Execute the following command on the HPC server:"
    cat <<'CMD'
cd ~/circRNA_agent
nohup snakemake \
    --snakefile workflow/Snakefile \
    --configfile config.yaml \
    --cores 36 \
    --resources mem_gb=300 \
    --use-singularity \
    --singularity-args "--bind /home3/choukaihsuan:/home3/choukaihsuan --bind /home3/choukaihsuan/reference:/reference" \
    --keep-going \
    --rerun-incomplete \
    > logs/pipeline_singularity.log 2>&1 &
echo "PID: $!"
CMD
    echo ""
    echo "NOTE: Set use_container: true in config.yaml before running."
    exit 0
fi

# ── Local (WSL2 / Linux) steps ────────────────────────────────────────────────

# Step 1: Build Docker image
echo "=== Step 1: Build Docker image ==="
# Must be run from the repository root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
cd "${REPO_ROOT}"

docker build -t "${IMAGE_NAME}:${IMAGE_TAG}" .
echo "Build complete: ${IMAGE_NAME}:${IMAGE_TAG}"

# Step 2: Local smoke test
echo ""
echo "=== Step 2: Local smoke test ==="
docker run --rm "${IMAGE_NAME}:${IMAGE_TAG}" CIRIquant --version
docker run --rm "${IMAGE_NAME}:${IMAGE_TAG}" snakemake --version
docker run --rm "${IMAGE_NAME}:${IMAGE_TAG}" Rscript -e "library(edgeR); cat('edgeR OK\n')"
echo "Smoke test passed."

# Step 3: Push to Docker Hub
echo ""
echo "=== Step 3: Push to Docker Hub ==="
echo "Make sure you are logged in:  docker login"
docker tag "${IMAGE_NAME}:${IMAGE_TAG}" "${FULL_IMAGE}"
docker push "${FULL_IMAGE}"
echo "Pushed: ${FULL_IMAGE}"

echo ""
echo "====================================================="
echo " Done. Next steps on the HPC server:"
echo "   bash containers/build_and_deploy.sh --hpc"
echo " Then set use_container: true in config.yaml"
echo "====================================================="
