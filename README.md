<h2 align="center">Boolean Gene</h2>
A Python package for analyzing Boolean gene networks.

## What are Boolean gene networks?
Interactions between genes and their regulators (such as transcription factors or RNA-binding proteins) provide the machinery for cells to respond to environmental stimuli by adopting an appropriate transcriptional response. The set of interconnected gene-regulator interactions – known as a gene regulatory network (GRN) – can be modeled to predict the cellular response to a signal. However, GRNs arising in nature are highly complex systems with dozens to hundreds of genetic regulators and effectors, making their full characterization computationally intractable. To nonetheless study their behavior, GRNs can be modeled as Boolean networks: each gene is represented by a node whose on/off state is updated at every time step according to a Boolean function of the gene’s regulators. A straightforward example is Elowitz and Leibler’s repressilator (1), where three genes are connected by NOT gates in a unidirectional loop. This GRN is designed to oscillate, and the Boolean model captures this behavior. For instance, when loaded with an initial condition where only one gene is on, the “on” state rotates unidirectionally through the loop (Fig. 1). 

<img src="/assets/repressilator1.png" width="30%"></img> <img src="/assets/repressilator2.png" width="30%"></img> <img src="/assets/repressilator3.png" width="30%"></img>  
**Figure 1:** Boolean model of the repressilator captures its oscillatory behavior. Genes are represented by nodes, and repressing interactions are represented by NOT gates.

## What does this software do?
This package allows users to perform stability analysis and measure topological features (e.g., loops or clusters) of an input Boolean gene network. The input can be either a user-defined network or a pre-defined network from the CellCollective database (2). A strength of this package is its ability to generate randomized networks that nevertheless preserve the logical structure of the input network, including the in-degrees, out-degrees, and Boolean function of each gene. This feature is integrated into the stability analysis pipeline such that users can compare the input network to large numbers of random networks that are logically identical but topologically different from the input. This analysis allows for the characterization of topological modifiers of GRN stability, which can have applications in synthetic biology, systems biology, and biomedical drug discovery. 

## Installation
```
pip install booleangene
```

## Usage
This package contains functions for running stability analysis, counting a network's loops, and measuring the clustering coefficient of each gene in a network (use `booleangene -h` to view all available functions). Each of these functions generates a user-defined number of randomized networks that preserve the Boolean functions of the original input network, allowing for comparison to a null model.

For instance, to run the stability analysis function:
```
booleangene run_stability_analysis -i "/path/to/input.txt" -o "/path/to/output/directory" 
```  
where the input .txt file should define each gene's associated Boolean function using AND, OR, and NOT connectives. For example inputs, see the text files in the `sample_inputs` directory, which were downloaded from the Cell Collective database (2).

## Citations
1. Elowitz, M., Leibler, S. A synthetic oscillatory network of transcriptional regulators. Nature 403, 335–338 (2000).
2. Helikar, T., Kowal, B., McClenathan, S. et al. The Cell Collective: Toward an open and collaborative approach to systems biology. BMC Syst Biol 6, 96 (2012).
