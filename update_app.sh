#!/bin/bash
set -euo pipefail
#Change to the scGenesElite directory
cd ~/scGenesElite
git fetch origin Master

# Check if the local copy of the application is up to date
if ! git diff --quiet HEAD origin/Master; then
  # If not, pull the latest version from Github
  git pull --ff-only origin Master
  
  # Create a fresh Docker image
  docker build -t scgenes .

  # Stop the old Docker container
  docker stop elite

  # Start a new Docker container
  docker run -d --rm --name elite -p 8447:3838 scgenes
  
  # remove the old image
  docker image prune -f
fi
