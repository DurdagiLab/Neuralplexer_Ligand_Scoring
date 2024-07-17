# Neuralplexer_Ligand_Scoring

This script generates a pipeline for finding the rank 1 ligand and the related protein structure generated by the NeuralPlexer generative model. It scores the ligand binding affinity using Schrodinger Maestro.

## Requirements

- **NeuralPlexer Generative Model**
  - Install from their [GitHub repository](https://github.com/zrqiao/NeuralPLexer).
  - For more information, refer to their article: [DOI: 10.1038/s42256-024-00792-z](https://doi.org/10.1038/s42256-024-00792-z).
- **Schrodinger Maestro**
  - Used for the scoring method. (This script uses Schrodinger version 2018)
- **Python Packages**
  - numpy: `pip install numpy`
  - rdkit: `pip install rdkit`
  - shutil: `pip install shutil`
  - glob: `pip install glob`

## Usage

1. Install all the required software and packages listed above.
2. Run NeuralPlexer to generate the ligand and protein structure outputs.
3. Update the script to change the path to the NeuralPlexer output directory.
