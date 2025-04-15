CONDA_BASE_DIR="/root/miniconda"
CONDA_SH_PATH="$CONDA_BASE_DIR/etc/profile.d/conda.sh"

# Lazy load conda to improve shell startup time
conda() {
    # Remove the function after first execution
    unset -f conda
    
    # Initialize conda - simple version
    . "$CONDA_SH_PATH"
    
    # Execute the original command
    conda "$@"
} 