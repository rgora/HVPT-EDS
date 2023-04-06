# HVP-EDS helper scripts

xyz2eds.py
==========

xyz2eds.py is an open source script, written in Python, for preparing a set of input files for hybrid variational-perturbational energy decomposition scheme (HVP-EDS) implemented in my custom version of GAMESS(US).

geds.py
==========
geds.py is an open source script, written in Python, for analyzing a set of output files of HVP-EDS calculations.

## Support ##

Please contact the author [mailto:robert.gora@pwr.edu.pl Robert Góra].

## Tutorials & References ##

For the details of the approach and the present implementation I can recommend our recent papers [1,2].

The usage of the script itself is rather self-explanatory, however, I have prepared a short tutorial which is available through project's wiki.

## References ##

The hybrid variational–perturbational interaction energy decomposition scheme (VP-EDS) was proposed independently by Gutowski et al. [^1] and Sokalski et al. [^2] in the late 1980s. Historically the first quantum chemical theory of intermolecular interactions which dates back to 1930s was due to London et al. and was based on the stan- dard Rayleigh–Schrodinger perturbation expansion [3]. Even though its modern reformulation, namely the symmetry–adapted perturbation theory (SAPT) is nowadays probably the most advanced, elegant and successful approach in this field [3, 4], we believe that the VP–EDS scheme outlined below may still provide in certain cases a plausible alternative. Particularly in the case of extended systems, since due to its present implementation and interfacing with the GAMESS software package [5] and taking advantage of the direct SCF approach and multi-core architecture, it allows to perform calculations for systems which are often far beyond the reach of other similar approaches. 

### Partitioning of the interaction energy of a dimer ###

The total interaction energy of a dimer, calculated by a supermolecular approach using the second–order Møller–Plesset perturbation theory (MP2), is partitioned into the Hartree–Fock (HF) and the electron correlation interaction energy components:
'''math
    \Delta E^{\rm MP2} = \Delta E^{\rm HF} + \Delta E_{\rm corr}^{\rm MP2}.
'''
The HF term can be partitioned into the Heitler--London (HL) and the
delocalization ($\Delta E^{\rm HF}_{\rm del}$) energy components.  The former
involves the electrostatic interactions of unperturbed monomer charge
densities, $\epsilon_{\rm el}^{\rm (10)}$, as well as the associated exchange
repulsion $\Delta E_{\rm ex}^{\rm HL}$ and can be calculated after L{\"o}wdin
\citep{lovdin-ap-v5-1956-p1} as a difference between the expectation value of
the Hamiltonian, defined as a sum of the free--monomer Hamiltonians
$\mathcal{H}_i$ and $\mathcal{H}_j$ and the intermolecular perturbation
$\mathcal{V}_{ij}$, for the wavefuncion given as an antysymmetrized Hartree
product of monomers' wavefunctions and the sum of monomer energies:
\begin{equation}
\begin{split}
    \Delta E^{\rm HL} & = N_{ij}^{HL}\left\langle
    \mathcal{A}\Psi_i\Psi_j\right\rvert
    \mathcal{H}_i + \mathcal{H}_j + \mathcal{V}_{ij} \left\lvert
    \mathcal{A}\Psi_i\Psi_j\right\rangle \\
    & - \left( E_i^{HF} + E_j^{HF} \right)
\end{split}
\end{equation} 
where $N_{ij}^{HL}$ is the normalization constant.
The delocalization term is estimated as a difference between the HF and the HL
interaction energies and encompasses the induction and the associated exchange
effects due to the Pauli exclusion principle. 
\begin{equation}
\begin{split}
    \Delta E^{\rm HF} & = \Delta E^{\rm HL} + \Delta E^{\rm HF}_{\rm del} \\ 
                      & = \epsilon_{\rm el}^{\rm (10)} + \Delta E_{\rm ex}^{\rm HL}
                        + \Delta E^{\rm HF}_{\rm del}
\end{split}
\end{equation}
The second order electron correlation term, $\Delta E_{\rm corr}^{\rm MP2}$,
includes the second order dispersion interaction, $\epsilon_{\rm disp}^{\rm
(20)}$, as well as the electron correlation correction to the first order
electrostatic interaction, $\epsilon_{\rm el,r}^{\rm (12)}$, and the remaining
electron correlation effects ($\Delta E_{\rm ex}^{\rm (2)}$). The latter term
accounts mainly for the uncorrelated exchange--dispersion and electron
correlation corrections to the Hartree--Fock exchange
repulsion \citep{chaasiski-mp-v63-1988-p205, cybulski-jcp-v92-1990-p4357}.
\begin{equation}
    \Delta E_{\rm corr}^{\rm MP2} = \epsilon_{\rm el,r}^{\rm (12)} 
    + \epsilon_{\rm disp}^{\rm (20)} + \Delta E_{\rm ex}^{\rm (2)}
\end{equation} 
The $\epsilon_{\rm el}^{\rm (10)}$ and the $\epsilon_{\rm disp}^{\rm (20)}$
terms are obtained in the standard polarization perturbation
theory, whereas the $\epsilon_{\rm el,r}^{\rm (12)}$
term is calculated using the formula proposed by Moszy{\'n}ski {\it et
al.} \citep{moszyski-cpl-v166-1990-p609}. The indices in parenthesis denote perturbation
orders in intermolecular interaction operator and intramonomer correlation
operator respectively.
In order to account for higher--order electron correlation effects, this scheme
can be extended to higher orders of perturbation theory. However, it is
often augmented instead with the supermolecular electron correlation
corrections estimated using the Coupled--Clusters (CC) approach with singles,
doubles, and noniterative triples i.e. $\Delta E_{\rm corr}^{\rm CCSD}$
and $\Delta E_{\rm corr}^{\rm CCSD(T)}$.

### Partitioning of the interaction energy of many body complexes ###

The above scheme can be generalized also to many body systems 
\citep{chalasinski-cr-v94-1994-p1723,gora-jcp-v117-p1031,
gora-jpcb-v109-2005-p2027}. In the supermolecular approach the total
interaction energy of an $n$--fragment system can be partitioned into
$2,3,\ldots,n$ body components as follows:
%
\begin{equation}\label{eq:pn}
\begin{split}
    & \Delta E(12\ldots n) = E(12\ldots n) - \sum_{i=1}^{n} E(i) = \\
    & \quad = \sum_{i=1}^{n-1}\sum_{j>i}^{n} \Delta^2 E(ij) 
            + \sum_{i=1}^{n-2}\sum_{j>i}^{n-1}\sum_{k>j}^{n} \Delta^3 E(ijk) \\
    & \qquad + \ldots + \Delta^n E(12\ldots n).
\end{split}
\end{equation} 
In the above equation $E(i)$, $E(ij)$, $E(ijk)$ denote a given property of
monomer $i$, dimer $ij$ and trimer $ijk$, respectively. The individual 2- and
3--body contributions are given by:
\begin{equation}
\begin{split}
    & \Delta^2 E(12)  = E(12) - \sum_{i=1}^2 E(i) \\
    & \Delta^3 E(123) = E(123) - \sum_{i=1}^3 E(i)
                      - \sum_{i=1}^2\sum_{j>i}^3 \Delta^2 E(ij).
\end{split}
\end{equation} 
Taking into account that both the electrostatic components and the second
order dispersion energy are pairwise additive the total MP2 interaction
energy, of a complex can be partitioned as follows:
\begin{equation}
\begin{split}
    \Delta E^{\rm MP2} & = \sum_{i=2}^{n}\Delta^i E^{\rm HF} +
                           \sum_{i=2}^{n}\Delta^i E^{\rm MP2}_{\rm corr} = \\
                       & = \epsilon_{\rm el}^{(10)} +
                           \sum_{i=2}^{n}\Delta^i E_{\rm ex}^{\rm HL} +
                           \sum_{i=2}^{n}\Delta^i E_{\rm del}^{\rm HF} + \\
                       & + \epsilon_{\rm el,r}^{(12)} +
                           \epsilon_{\rm disp}^{(20)} +
                           \sum_{i=2}^{n}\Delta^i E_{\rm ex}^{(2)} 
\end{split}
\end{equation} 
In all necessary calculations the complex centered basis set is consistently
used and, therefore, the results are basis set superposition error free due to
the full counterpoise correction \citep{boys-mol_phys-v19-p553}. In the case of
many body complexes this procedure is equivalent to the site--site function
counterpoise method.

## References ##

[^1]: R. W. Góra, W. Bartkowiak, S. Roszak and J. Leszczyński, A New Theoretical Insight into the Nature of Intermolecular Interactions in the Molecular Crystal of Urea, J. Chem. Phys., 2002, 117, 1031–1039.
2 R. W. Góra, W. Bartkowiak, S. Roszak and J. Leszczyński, Intermolecular Interactions in Solution: Elucidating the Influence of the Solvent, J. Chem. Phys., 2004, 120, 2802–2813.
