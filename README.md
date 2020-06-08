# CS_foldrec

Requires PERL (v.5.8 or higher)

- FOLD RECOGNITION :: COMPLEMENTARITY SCORE
- CS_gl, CS_cp 
- Computed from side-chain shape complementarity (Sm) & electrostatic complementarity (Em) of buried amino acids at the globular protein interior 
- Burial Cutoff: 0.30
- Sequence -> Threaded onto a fold -> (Sm, Em) -> CSgl, CScp
- ---------------------------------------------------------------
- Reference: 
- ARTICLE| VOLUME 102, ISSUE 11, P2605-2614, JUNE 06, 2012
- Self-Complementarity within Proteins: Bridging the Gap between Binding and Folding
- Sankar Basu, Dhananjay Bhattacharyya, Rahul Banerjee
- Open ArchiveDOI:https://doi.org/10.1016/j.bpj.2012.04.029
- --------------------------------------------------------------

```Usage: ./CS_foldrec.pl <\$basename>.Sm <\$basename>.Em```

