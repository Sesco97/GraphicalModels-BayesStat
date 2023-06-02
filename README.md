# GraphicalModels-BayesStat

## Table of Contents 
- [Introduction](#introduction "Goto introduction")
- [Technologies](#technologies "Goto technologies")
- [Results](#results "Goto results")

## Introduction
The aim of this project is to develop Bayesian methods for the analysis of multivariate categorical data.
In particular, we are interested in inferring dependence relations between categorical variables, also accounting for possible heterogeneity related to latent clustering structures in the data.

We adopt graphical models to represent dependence relations between variables: specifically, a *graphical model* is a probabilistic model for a collection of random variables based on a graph structure.
A graph $\mathcal{G} = (V,E)$ is made up of a set of *nodes* $V$ (representing variables) and a set of *edges* $E$ (representing dependence relations between nodes).
Therefore it can be used to model conditional dependence structures between random variables.

We tackle the problem of identifying dependence relations between variables as a model (graph) selection problem.
Target of the analysis is to approximate a posterior distribution over the space of all possible graph structures given the data matrix.

Two different conjugate models to deal with categorical data will be presented, based on:

- **Multivariate categorical (multinomial) distributions** with Hyper Dirichlet prior;
- **Latent multivariate Gaussian distributions** with Hyper-Inverse-Wishart prior (inference via data augmentation).

We will show that in both cases the posterior is known up to a normalizing constant and we will present the implementation of a Metropolis-Hastings algorithm in order to make inference on $\mathcal{G}$ given the data.

**A complete report of the project is available [here](Report.pdf)**.

## Technologies
This project has been entirely created with:

- R 4.3.0.

## Results
An assessment of the performance of the two models for inference on the graph, namely the Multinomial Dirichlet and the Latent Normal Inverse-Wishart, has been carried out both on a simulated and real dataset.

### Simulation
We first created 20 categorical datasets starting from as many randomly generated decomposable graphs (with 6 nodes) and then we ran both algorithms on each of these datasets.
After that, we computed the *Structural Hamming Distance* between the original graph generating the data and the graph estimated from the resulting chain.
In particular, two ways of estimating the graph have been considered:

1. ***Median Probability Graph (MPG)***: the graph obtained by including all the edges that were in at least 50% of the graphs visited by the chain;
2. ***Maximum a Posteriori (MAP)***: the most recurrent graph in the chain.

The obtained results are shown in section 5.1 of the [report](Report.pdf).

### Real Dataset
Inference has been performed on the [*Congressional Voting Records*](https://archive.ics.uci.edu/ml/datasets/congressional+voting+records) dataset.
Such dataset includes votes for each of the U.S. House of Representatives Congressmen on the 16 key votes identified by the CQA.

In particular, we divided the congressmen into two groups (republicans and democrats) and we ran the Multinomial-Dirichlet model individually on each of them. Note that we chose the Multinomial-Dirichlet model mainly because of the lower computational time.

The obtained results are shown in section 5.2 of the [report](Report.pdf).
