#!/bin/bash
git clone https://hub.fastgit.org/Gaius-Augustus/Augustus.git
cd Augustus
docker build -t augustus:latest .