# ADMETlab models for evaluation of drug candidates

A series of models for the systematic ADMET evaluation of drug candidate molecules. Models include blood-brain barrier penetration; inhibition and substrate affinity for CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, and pgp; F 20% and F 30% bioavailability; human intestinal absorption; Ames mutagenicity; skin sensitization; plasma protein binding; volume distribution; LD50 of acute toxicity; human hepatotoxicity; hERG blocking; clearance; half-life; Papp (caco-2 permeability); LogD distribution coefficient; and LogS solubility.

## Identifiers

* EOS model ID: `eos2re5`
* Slug: `admetlab`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Classification, Regression`
* Output: `Experimental value`
* Output Type: `Float`
* Output Shape: `List`
* Interpretation: Regression models provide a numerical result (LogS (log mol/L), LogP (distribution coefficient), Papp (Caco-2 permeability in cm/s), PPB (%)). Classifications provide the probability of activity according to ADMETlab thresholds.

## References

* [Publication](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0283-x)
* [Source Code](https://github.com/ifyoungnet/ADMETlab)
* Ersilia contributor: [svolk19-stanford ](https://github.com/svolk19-stanford )

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos2re5)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2re5.zip)

## Citation

If you use this model, please cite the [original authors](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0283-x) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a GPL-3.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!