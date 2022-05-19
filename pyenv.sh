#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2022-05-19 13:55:49
# @DESCRIPTION:

# Number of input parameters



[ -d venv ] || eval '`which python3` -m venv venv'

source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt