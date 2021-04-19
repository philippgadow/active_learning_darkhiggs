#!/bin/bash

if [ -d "venv/" ]; then
  ### Take action if $DIR exists ###
  echo "Using python virtual environment in venv/..."
else
  ###  Control will jump here if $DIR does NOT exists ###
  echo "Info: venv/ not found. Setting up python virtual environment for the first time..."
  wget https://files.pythonhosted.org/packages/33/bc/fa0b5347139cd9564f0d44ebd2b147ac97c36b2403943dbee8a25fd74012/virtualenv-16.0.0.tar.gz
  tar xvf virtualenv-16.0.0.tar.gz
  python virtualenv-16.0.0/virtualenv.py venv
  rm -r virtualenv-16.0.0  virtualenv-16.0.0.tar.gz
  pip install -r requirements.txt
  echo "Python virtual environment set up! Activating it now..."
fi

source venv/bin/activate
