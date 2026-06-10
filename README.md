# A quantitative portrait of habituation in *Stentor coeruleus*

**Tejas Ramdas**<sup>1*</sup>, **Nhi Doan**<sup>2</sup>, **Austen Theroux**<sup>2</sup>, **Samuel J. Gershman**<sup>2,3</sup>

<sup>1</sup> Program in Neuroscience, Harvard University  
<sup>2</sup> Department of Psychology, Harvard University  
<sup>3</sup> Center for Brain Science, Harvard University  
<sup>*</sup> Correspondence: tejasramdas@g.harvard.edu

## Abstract

Habituation—the decrement in response to a series of stimuli—is a widespread form of learning observed across many organisms, including the unicellular organism *Stentor coeruleus*. A lesser-known feature of *Stentor* habituation, shared with animals, is potentiation: faster habituation to a second stimulus series despite partial or complete recovery of responsiveness before that series begins. This suggests that although the first-order habituation memory can decay during the recovery period between the two series, a persistent second-order memory mediates faster relearning. We investigate the response profile of *Stentor* across a range of stimulation frequencies and recovery periods to identify the timescales at which these memory traces operate. We introduce a statistical framework to infer both population and single-cell learning parameters, allowing us to quantify prior qualitative findings and examine relationships among parameters across cells. Two key findings are that potentiation is frequency-sensitive, and that recovery and potentiation are decoupled, consistent with a serial and hierarchical cascade of leaky integrator units underlying these processes. This quantitative portrait provides a foundation for mechanistic modeling of intracellular memory in *Stentor*.

## Code

This repo contains code and data files for reproducing the analysis and figures in the paper.

- `paper.jl` — main script for generating all figures
- `processed_data/collated.jld2` — collated behavioral data
- `processed_data/inferred_chains.jld2` — pre-computed MCMC chains

## Citation

Ramdas, T., Doan, N., Theroux, A., & Gershman, S. J. (2026). A quantitative portrait of habituation in *Stentor coeruleus*. bioRxiv. https://doi.org/10.64898/2026.06.09.731162

```bibtex
@article {Ramdas2026,
	author = {Ramdas, Tejas and Doan, Nhi and Theroux, Austen and Gershman, Samuel J},
	title = {A quantitative portrait of habituation in Stentor coeruleus},
	elocation-id = {2026.06.09.731162},
	year = {2026},
	doi = {10.64898/2026.06.09.731162},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2026/06/10/2026.06.09.731162},
	eprint = {https://www.biorxiv.org/content/early/2026/06/10/2026.06.09.731162.full.pdf},
	journal = {bioRxiv}
}
```
