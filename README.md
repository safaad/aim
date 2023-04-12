# AIM: Alignment-in-Memory A Framework for High-throughput Sequence Alignment using Real Processing-in-Memory Systems
AIM is a framework for high-throughput pairwise sequence alignment
using processing-in-memory. It targets [UPMEM](https://www.upmem.com/), the first publicly-available
general-purpose programmable PIM architecture. AIM dispatches a large number
of sequence pairs across different memory modules and aligns each pair
using compute cores within the memory module where it resides. It supports multiple alignment algorithms including NW, SWG, GenASM [[1](#myfootnote1)], WFA [[2](#myfootnote2)], and WFA-
adaptive, and includes two implementations of each algorithm that
manage the UPMEM memory hierarchy differently and are suitable
for different read lengths.

We evaluate AIM on a real UPMEM system and compare the throughput it can achieve with that achieved by server-grade multi-threaded CPU systems running at full scale.
Our evaluation shows that a real PIM system can substantially outperform CPU systems for a wide variety of algorithms, read lengths, and edit distance thresholds. For example, for WFA-adaptive, the state-of-the-art sequence alignment algorithm, AIM achieves a speedup of up to 2.56x when data transfer time is included, and up to 28.14x when data transfer time is not included.

## Citation
For more information on this project you can refer to the following papers:

Long paper:
> Safaa Diab, Amir Nassereldine, Mohammed Alser, Juan Gómez Luna, Onur Mutlu, Izzat El Hajj, "**[A Framework for High-throughput Sequence Alignment using Real Processing-in-Memory Systems"](https://doi.org/10.1093/bioinformatics/btad155)**, Bioinformatics, 2023; btad155, https://doi.org/10.1093/bioinformatics/btad155


HiCOMB22 short paper ([slides](https://people.inf.ethz.ch/omutlu/pub/WFA-PairwiseAlignment-in-PIM_hicomb22-GPU-hicomb22-talk)):

> Safaa Diab, Amir Nassereldine, Mohammed Alser, Juan Gomez Luna, Onur Mutlu, & Izzat El Hajj (2022). "**[High-throughput Pairwise Alignment with the Wavefront Algorithm using Processing-in-Memory](https://www.computer.org/csdl/proceedings-article/ipdpsw/2022/974700a163/1Fu98na0V3y)**". In 2022 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW). IEEE.

If you find this project useful, please cite these papers.
## Repository Structure
Each folder of this repository has the PIM implementations relative to each alignment algorithm. For each alignment algorithm, we provide two PIM implementations: (1) DPU-WRAM implementation stores the alignment data in the WRAM, and (2) DPU-MRAM implementation stores the alignment data in the MRAM and uses the WRAM as cache memory.

GenASM's PIM implementations are found in this submodule of AIM framework https://github.com/safaad/aim-genasm
```bash
├───Datasets
├───NW
│   ├───DPU-MRAM
│   └───DPU-WRAM
├───SWG
│   ├───DPU-MRAM
│   └───DPU-WRAM
└───WFA
    ├───DPU-MRAM
    └───DPU-WRAM
```
## Instructions

### Prerequisites
AIM is designed to run on a server with a real UPMEM-PIM modules. However, the functional simulator included in the [UPMEM SDK](https://sdk.upmem.com/) can be also used to run the implementations. The used SDK version in AIM is [2021.3.0](https://sdk.upmem.com/).

### Getting Started
Clone the repository:
```bash
git clone --recurse-submodules https://github.com/safaad/aim.git

# Or if you don't want to include aim-genasm
git clone https://github.com/safaad/aim
cd aim
```
### Running AIM
AIM provides a python script to run each PIM implementation `run_*_pim.py`. The script estimates the upper limit of the used WRAM and MRAM memory in order to determine the number of DPU threads needed `NR_TASKLETS`.It takes as an input features of the read pairs dataset such as read length, edit distance threshold, and the number of sequence pairs to align. Below is an example of running WFA MRAM implementation (same method can be used for all versions) while enabling backtracing on a dataset with sequence length 100bp, ED=1%, and 40K sequence pairs:
```bash
cd WFA/DPU-MRAM

# Script usage
> python run-wfa-pim-mram.py
> usage: run-wfa-pim-mram.py [-h] -i INPUT [-o OUTPUT] -l READ_LENGTH -e ERROR
                           -n NUMBER_READS [-m MATCH_COST] [-x MISMATCH_COST]
                           [-g GAP_OPENING] [-a GAP_EXTENDING] [-b] [-r]
                           [-t NR_OF_TASKLETS] [-d NR_OF_DPUS]
                           
# Example of running WFA MRAM implementation
python run-wfa-pim-mram.py -i ../../Datasets/sample-l100-e1-40K.01 -l 100 -e 0.01 -n 40000 -b -d 2500

```
Note that this script uses heuristics to estimate the upper limit of the used WRAM and MRAM memory, so it might underestimate the `NR_TASKLETS` to avoid running out of memory.


Or, the user can compile and run these implementations manually but needs to determine the `NR_TASKLETS` and the `WRAM_SEGMENT` size assigned to each thread if the implementation uses the custom dynamic memory allocator (WFA and GenASM). Also, the sequence length, max alignment score, and the number of sequence pairs should be added at the compile and run time. For example:
```bash
cd WFA/DPU-MRAM

# Compile example
make NR_TASKLETS=19 NR_DPUS=2500 FLAGS="-DMAX_SCORE=25 -DREAD_SIZE=112 -DBACKTRACE -DWRAM_SEGMENT=2122"

# Run Example
./build/host ../../Datasets/sample-l100-e1-40K.01 ./out 40000
```
Each line of the output file will contain the number of the aligned read-reference pair, the alignment score (edit distance in case of GenASM), and the CIGAR string if the backtracing is enabled.

## Contact

For further questions and suggestions, feel free to reach out syd04@aub.edu.lb

## References
* <a name="myfootnote1">[1] </a> Damla Senol Cali, Gurpreet S. Kalsi, Zülal Bingöl, Can Firtina, Lavanya Subramanian, Jeremie S. Kim, Rachata Ausavarungnirun, Mohammed Alser, Juan Gomez-Luna, Amirali Boroumand, Anant Nori, Allison Scibisz, Sreenivas Subramoney, Can Alkan, Saugata Ghose, and Onur Mutlu.
["GenASM: A High-Performance, Low-Power Approximate String Matching Acceleration Framework for Genome Sequence Analysis."](https://people.inf.ethz.ch/omutlu/pub/GenASM-approximate-string-matching-framework-for-genome-analysis_micro20.pdf)
In _Proceedings of the 53rd International Symposium on Microarchitecture (MICRO),_ Virtual, October 2020.


* <a name="myfootnote2">[2] </a> Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.




