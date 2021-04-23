FROM atlas/analysisbase:latest
ADD . /opt/SimpleAnalysis
# set a bunch of environment variables based on CI
ARG CI_COMMIT_SHA
ARG CI_COMMIT_REF_SLUG
ARG CI_COMMIT_TAG
ARG CI_JOB_URL
ARG CI_PROJECT_URL
ENV CI_COMMIT_SHA=$CI_COMMIT_SHA
ENV CI_COMMIT_REF_SLUG=$CI_COMMIT_REF_SLUG
ENV CI_COMMIT_TAG=$CI_COMMIT_TAG
ENV CI_JOB_URL=$CI_JOB_URL
ENV CI_PROJECT_URL=$CI_PROJECT_URL
# actually build, install, and add to release_setup.sh
RUN source /release_setup.sh && \
    cmake -S /opt/SimpleAnalysis -B /tmp/SimpleAnalysis -DCMAKE_INSTALL_PREFIX=/usr/SimpleAnalysis && \
    cmake --build /tmp/SimpleAnalysis -j4 && \
    sudo /usr/local/bin/cmake --install /tmp/SimpleAnalysis && \
    rm -rf /tmp/SimpleAnalysis && \
    sudo chmod 666 /release_setup.sh && \
    sudo printf '\n# Set up the SimpleAnalysis code\nsource /usr/SimpleAnalysis/setup.sh\necho "Configured SimpleAnalysis from: $SimpleAnalysisProject_DIR"' >> /release_setup.sh && \
    sudo chmod 644 /release_setup.sh
