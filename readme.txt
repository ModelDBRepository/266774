March 18, 2021: Initial code upload (Markovian FF and Non-Markovian)

January 24, 2023: -Changes made to Markovian code: (upload of Markovian A2A)
		*Original Markovian code (Markovian FF) incorrectly restricted Messenger to Timer intercolumnar 
		connections were such that they were limited to be feed-forward (i.e. ordinal to the presented sequence)
		*Updated Markovian code (Markovian A2A) correctly allows for any (all-to-all or A2A) Messenger to Timer intercolumnar
		connections. To allow for this, the following changes were made:
			*Feed-forward threshold changed from 20Hz to 30Hz
			*T^d,ff_max changed from .00345 to .0045 
			*eta^p_ff changed from 20 x 3500 ms^-1 to 8.8 x 3500 ms^-1
			*eta^d_ff changed from 15 x 3500 ms^-1 to 10 x 3500 ms^-1
			*eta_ff changed from 0.25 to 0.4
		*Supplementary Table 1 has been updated include values for both Markovian FF and Markovian A2A
	      -Updates made to Supplementary Tables:
		*Original Supplementary Tables contained incorrectly scaled or incorrectly stated parameters.
		Please see[1] for more details. The Supplementary Tables included with this upload have been 
		updated to rectify these errors.

The above errors were found during a replication study (Zajzon et al. 2023[1]), which also proposed the solution of 
raising the feed-forward threshold (implemented in Markovian A2A). We thank the authors of that study for their contributions to these changes.

Additionally, Zajzon et al. 2023[1] produced a well-documented and effecient implementation of this model in NEST, so we encourage those 
interested in using the model for their own purposes to also consider their implementation as a potential code base[2].


1.Zajzon, B., Duarte, R. & Morrison, A. Towards reproducible models of sequence learning: replication and analysis of a modular spiking network with reward-based learning. bioRxiv 2023.01.18.524604 (2023) doi:10.1101/2023.01.18.524604.
2.Zajzon, B. Towards reproducible models of sequence learning: replication and analysis of a modular spiking network with reward-based learning. https://github.com/zbarni/re_modular_seqlearn (2022).

