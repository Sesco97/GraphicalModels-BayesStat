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

**A complete report of the project is available [here]()**

## Technologies

## Results
