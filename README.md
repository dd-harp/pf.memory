# `pf.memory` <br><br> The Age of Infection, Immunity & Malaria Epidemiology 

This software is a set of computational algorithms that implements a synthetic approach to malaria epidemiology, defined narrowly to include exposure to the bites of infective mosquitoes, the dynamics of infection (including superinfection).  

The epidemiology of *Plasmodium falciparum* malaria presents a unique set of challenges due to the complex dynamics of infection, immunity, disease and infectiousness as well as treatment and chemo-protection, diagnostics and detection. Malaria can be measured in a dozen different ways, but it has been difficult to present a simple synthesis of malaria infection and disease in terms of the metrics that are commonly used in research and clinical surveillance. An important metric is the *Plasmodium falciparum* parasite rate, or *Pf*PR, defined as the average prevalence of malaria taken from a cross-sectional survey. Another metric, often measured as a covariate in research studies, is a parasite count, the number of parasites in a blood slide field counted by a light microscopist. In an old data set, collected during malariotherapy, parasite counts fluctuated substantially over the time course of an infection, but they were strongly statistically correlated with the *age of infection* or AoI (Henry JM, *et al.*, 2022)^[Henry JM, Carter A & Smith DL (2022) Infection age as a predictor of epidemiological metrics for malaria. Malar J 21; 117, https://doi.org/10.1186/s12936-022-04134-5]. In malaria, the *Pf*PR in several old studies had a characterstic shape when plotted against age. The patterns also differ by diagnostic method, by season, by sex, and by location. Parasite densities have been used as a correlate of disease.   


With so many interacting factors, it was a challenge to develop a model that could deal with everything. One approach to studying malaria infection has been to develop mechanistic models for the dynamics of malaria infection within a single host. The most prominent models of this type are OpenMalaria and eMod, but there have been several others. These computational models made it possible to develope comprehensive individual-based simulation models, or IBMs, of malaria for policy. While these approaches have been able to replicate the patterns, the outputs of the models are usually just as complex as the data collected in field studies. A synthesis of malaria epidemiology has proven elusive. 

We present the computational algorithms that support a probabilistic approach to malaria epidemiology.  We start with a semi-Markovian model of malaria exposure and infection, whose states are represented by random variables that describe the multiplicity of infection (MoI) in a host and the age of infection (AoI). Assuming that parasite densities are, at least statistically, predicted by the AoI, we can compute probability distribution functions describing parasite densities, parasite counts, and detection in an individual chosen at random from the population. From this, we present a model for parasite detection and parasite counts. This same approach has been extended to predict disease, immunity, treatment with anti-malarial drugs, and a brief period of chemo-protection. 

The probabilistic approach is both highly realistic and descriptive, but our goal was a synthesis. This synthesis involves a few steps:

1. We derive variables that can track the mean MoI, the mean AoI and all its moments, and the probability of detection. 

2. Hybrid models for the MoI already exist, but we show how this approach can be extended to systems of differential equations that track the mean and higher order moments of the distribution of the AoI. 

3. We a new random variable describing the age of the youngest infection (AoY). We show how the variable serves as a basis for computing parasite density distributions in complex infections. 

4. We derive a hybrid model with a variable that approximates the mean AoY. 

5. We demonstrate that a simple system of ordinary differential equations can be used to in place of the probabilistic models to compute the values of the research metrics. 

To put it another way, we can reduce the behavior of these highly complex probabilistic systems to a simple system of equations that has a high degree of accuracy. The computational and conceptual simplicity of hybrid models have some simplicity over compartmental models and stochastic individual-based models, and with the supporting probabilistic framework, provide a sound basis for a synthesis of observational malaria epidemiology. 


| A | B  |
| --- | --- |
| a  | b  |


## References


