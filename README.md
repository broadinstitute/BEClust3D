# BEClust3D

## Colab Notebook Example

Examples of how to use the BE method development pipeline is provided in the following colab notebook: 

| Link | Description |
|---------|-------------|
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pia-francesca/ema/blob/main/colab_notebooks/application_example_ion_channel_proteins.ipynb) | Example of how to use the BE pipeline for single screen input

## 🚀 Installation Guide

### 🐧 Linux (x86_64)

#### ✅ Option 1: Using Conda Environment YAML
Create the environment directly from the pre-configured YAML file:

```bash
conda env create -f requirements_linux64.yaml

```
####  Option 2: Manual Setup with `requirements.txt`
Create an anaconda environment
``` bash
conda create -n beclust3d python=3.10
conda activate beclust3d
```
Install external dependencies
```bash
conda install salilab::dssp
conda install -c conda-forge boost=1.73 
conda install bioconda::clustalo
```
Install Python dependencies
```bash
pip install -r requirements.txt
```

### For Google Colab
Add the following line at the top of your notebook:
```
! apt-get update
! apt-get install dssp clustalo
```
