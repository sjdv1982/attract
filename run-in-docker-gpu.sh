set -u -e

# Runs a command inside the ATTRACT Docker image
# - Current directory is mounted to /cwd, and the command is executed there

docker run --rm \
  -v `pwd`:/cwd \
  --workdir /cwd \
  --user $(id -u):$(id -g) \
  --shm-size 4GB \
  --gpus all \
  -t -i \
  rpbs/attract $*
