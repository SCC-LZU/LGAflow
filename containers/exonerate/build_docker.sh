#!/bin/bash
git clone https://hub.fastgit.org/nathanweeks/exonerate.git
cd exonerate
docker build -t exonerate .