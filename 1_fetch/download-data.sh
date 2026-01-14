#!/usr/bin/env bash

uvx geofetch -i sample-table.txt --add-convert-modifier --discard-soft -m `pwd` \
    && export SRARAW=`pwd` \
    && mkdir -p fastqs \
    && export SRAFQ=`pwd`/fastqs \
    && uvx looper run --config ./sample-table/looper_config.yaml -p local --output-dir . \
    && echo "complete"
