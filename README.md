# circRNA_agent

## Container Deployment

The pipeline ships a Docker image that can be converted to a Singularity `.sif`
for portable execution on any HPC cluster.

### 1 · Build and test locally (WSL2 / Linux)

```bash
# Run from the repo root (Docker Hub user is already set to choukaihsuan):
bash containers/build_and_deploy.sh
```

This will:
- Build `circrna-pipeline:1.0.0` from the `Dockerfile`
- Run a smoke-test (CIRIquant, snakemake, R/edgeR)
- Tag and push to Docker Hub

### 2 · Convert to Singularity on the HPC server

```bash
# On 172.16.0.178 (CentOS 7):
bash ~/circRNA_agent/containers/build_and_deploy.sh --hpc
```

This pulls the Docker Hub image and saves it as
`~/circRNA_agent/containers/circrna-pipeline_1.0.0.sif`.

The `--bind` flags map HPC host paths into the container:

| Host path | Container path | Contents |
|-----------|---------------|----------|
| `/home3/choukaihsuan` | `/home3/choukaihsuan` | home, results, sra_cache |
| `/home3/choukaihsuan/reference` | `/reference` | hg19 genome + indices |

### 3 · Run the pipeline with Singularity

```bash
# Set use_container: true in config.yaml first, then:
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
```

### When to set `use_container: true`

Set this flag whenever running with `--use-singularity`.  It switches the
CIRIquant config to `config/ciriquant_container.yaml`, which uses tool names
(no absolute paths) and the container's Java at
`/usr/lib/jvm/default-java/bin/java`.
