# ADMETlab
## Model identifiers
- Slug: admetlab
- Ersilia ID: eos2re5
- Tags: ADMET, bioactivity, toxicity

# Model description
A series of models for the systematic ADMET evaluation of drug candidate molecules. Models include blood-brain barrier penetration; inhibition and substrate affinity for CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, and pgp; F 20% and F 30% bioavailability; human intestinal absorption; Ames mutagenicity; skin sensitization; plasma protein binding; volume distribution; LD50 of acute toxicity; human hepatotoxicity; hERG blocking; clearance; half-life; Papp (caco-2 permeability); LogD distribution coefficient; and LogS solubility. 

- Input: SMILES 
- Output: A vector of classification and regression outputs. Each classification model contributes two values to the vector: a binary label value <model name>_label (e.g. whether a molecule is a CYP1A2 inhibitor or not) and float probability value <model name>_prob (e.g. the probability of a molecule acting as a CYP1A2 inhibitor). All classification label thresholds are at p=50%. Each regression model contributes one prediction float value <model name>_pred to the vector. The numerical ranges and meanings of the outputs are described in detail on the [ADMETlab website](https://admet.scbdd.com/home/interpretation/) and in the [ADMETlab publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0283-x).
- Model type: Classification and regression
- Training set: Training datasets for each model are documented on the ADMETlab website [here](https://admet.scbdd.com/home/interpretation/). 
- Mode of training: Pretrained. Model checkpoints were downloaded from the [GitHub repo](https://github.com/ifyoungnet/ADMETlab) provided by ADMETlab.

# Source code
Jie Dong, Ning-Ning Wang, Zhi-Jiang Yao, Lin Zhang, Yan Cheng, Defang Ouyang, Ai-Ping Lu, Dong-Sheng Cao. ADMETlab: a platform for systematic ADMET evaluation based on a comprehensively collected ADMET database. Journal of Cheminformatics, 2018, 10:29
- Code: https://github.com/ifyoungnet/ADMETlab/tree/master/example
- Checkpoints: https://github.com/ifyoungnet/ADMETlab

The plasma protein binding, volume distribution, clearance, half-life, hERG blocking, LD50 of acute toxicity, Papp, logD, and logS models take chemical descriptors as inputs. These descriptors are calculated using the no longer maintained ChemoPy (pychem) library archived [here](https://code.google.com/archive/p/pychem/source) and originally downloaded for use from [this repository](https://github.com/salotz/chemopy). 

# License
The GPL-v3 license applies to all parts of the repository that are not externally maintained libraries. This repository uses the externally maintained library "ADMETlab", located at [/model/framework](https://github.com/ersilia-os/eos2re5/tree/main/model/framework/ADMETlab) and licensed under a [GNU Affero General Public License](https://github.com/ersilia-os/eos2re5/blob/main/model/framework/ADMETlab/LICENSE).

# History 
- Model checkpoints were downloaded from https://github.com/ifyoungnet/ADMETlab on 7/28/2022.
- No modifications were made to the original code, though the code in [/example/run.py](https://github.com/ersilia-os/eos2re5/tree/main/model/framework/ADMETlab/example/run.py) was copied and adapted for use in [/code/python_2_main.py](https://github.com/ersilia-os/eos2re5/blob/main/model/framework/code/python_2_main.py)
- The ChemoPy chemical descriptor library was downloaded from https://github.com/salotz/chemopy on 8/24/2022. Each import reference to the legacy pybel library in the ChemoPy source code was changed to instead import pybel from the openbabel library.

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission or [volunteer](https://www.ersilia.io/volunteer) with us!
