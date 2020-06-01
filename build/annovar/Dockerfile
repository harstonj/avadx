FROM kfdrc/annovar:latest as base

# setup bio-node
LABEL bio-node=v1.0

# set environment variables
RUN ln -s "$TOOLS/$TOOL_NAME/current/bin" /app
WORKDIR /app

# set app ENTRYPOINT
ENTRYPOINT ["/bin/bash"]

# set app CMD
CMD []
