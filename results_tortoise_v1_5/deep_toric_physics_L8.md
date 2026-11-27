# Theorem L8 — Deep Toric Memory Spectrum, Algebraic Null-Tower & Evolutionary Slowdown Clock (v1.5.0)

> **Provenance:** `2026-11-27T08:21:20Z` · **30 Species** · `zero_rng = true` · Zero free fitting parameters.

## 1. Multi-Lag Euchromatin Memory Decay (`d = 3..90 bp`), Null-Tower & Intra-Chromosome Split-Half Reliability

| Species | Life | I_nonbr(3) | I_nonbr(10) | I_nonbr(30) | I_nonbr(90) | R_mem(30/3) | Beyond-M1 % | Beyond-M2 % | Pal Share % | SplitHalf cos(M2) | SplitHalf Δ3600D |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Immortal_Jellyfish** | ∞ | 0.0060 | 0.0023 | 0.0012 | 0.0005 | 0.194 | 33.2% | 20.0% | 3.97% | 0.9900 | 0.17 pp |
| **Hydra_Immortal** | ∞ | 0.0332 | 0.0269 | 0.0190 | 0.0061 | 0.573 | 39.2% | 55.7% | 4.67% | 0.9990 | 0.17 pp |
| **Ocean_Quahog** | 507 yr | 0.0105 | 0.0060 | 0.0032 | 0.0013 | 0.304 | 63.4% | 45.8% | 3.29% | 0.9930 | 0.08 pp |
| **Greenland_Shark** | 400 yr | 0.0132 | 0.0081 | 0.0069 | 0.0048 | 0.525 | 48.9% | 36.4% | 2.59% | 0.9849 | 1.45 pp |
| **Rougheye_Rockfish** | 205 yr | 0.0208 | 0.0109 | 0.0047 | 0.0018 | 0.225 | 67.7% | 61.0% | 3.15% | 0.9852 | 0.88 pp |
| **Aldabra_Tortoise** | 190 yr | 0.0084 | 0.0053 | 0.0037 | 0.0024 | 0.443 | 38.2% | 26.1% | 3.05% | 0.9942 | 0.80 pp |
| **Galapagos_Tortoise** | 175 yr | 0.0081 | 0.0049 | 0.0034 | 0.0021 | 0.413 | 38.3% | 26.3% | 3.15% | 0.9946 | 0.92 pp |
| **Human_T2T** | 122.5 yr | 0.0262 | 0.0194 | 0.0134 | 0.0053 | 0.514 | 60.0% | 63.1% | 1.58% | 0.9944 | 2.07 pp |
| **Blue_Whale** | 110 yr | 0.0185 | 0.0122 | 0.0078 | 0.0049 | 0.421 | 50.0% | 58.6% | 2.23% | 0.9910 | 1.09 pp |
| **American_Lobster** | 100 yr | 0.0452 | 0.0299 | 0.0248 | 0.0110 | 0.549 | 78.3% | 121.8% | 1.96% | 0.9289 | 2.17 pp |
| **Desert_Tortoise** | 80 yr | 0.0087 | 0.0055 | 0.0039 | 0.0026 | 0.443 | 38.0% | 26.4% | 3.06% | 0.9941 | 1.12 pp |
| **American_Alligator** | 75 yr | 0.0108 | 0.0068 | 0.0046 | 0.0032 | 0.425 | 33.4% | 27.9% | 2.90% | 0.9970 | 0.31 pp |
| **African_Elephant** | 70 yr | 0.0133 | 0.0093 | 0.0050 | 0.0026 | 0.374 | 41.6% | 41.2% | 1.66% | 0.9809 | 3.48 pp |
| **Painted_Turtle** | 42 yr | 0.0099 | 0.0062 | 0.0044 | 0.0030 | 0.442 | 39.5% | 29.8% | 3.08% | 0.9957 | 0.77 pp |
| **Naked_Mole_Rat** | 37 yr | 0.0139 | 0.0089 | 0.0057 | 0.0032 | 0.406 | 46.1% | 46.3% | 1.92% | 0.9818 | 1.31 pp |
| **Greater_Mouse_Eared_Bat** | 37 yr | 0.0210 | 0.0138 | 0.0085 | 0.0038 | 0.403 | 48.7% | 66.7% | 3.24% | 0.9914 | 1.53 pp |
| **Black_Garden_Ant** | 28.8 yr | 0.0218 | 0.0138 | 0.0079 | 0.0038 | 0.362 | 55.1% | 58.4% | 5.03% | 0.9675 | 2.52 pp |
| **Softshell_Turtle** | 25 yr | 0.0134 | 0.0071 | 0.0047 | 0.0032 | 0.352 | 47.4% | 37.9% | 3.14% | 0.9945 | 1.77 pp |
| **Pacific_Oyster** | 20 yr | 0.0133 | 0.0075 | 0.0029 | 0.0007 | 0.220 | 52.1% | 40.2% | 4.61% | 0.9867 | 1.11 pp |
| **Dampwood_Termite** | 15 yr | 0.0100 | 0.0036 | 0.0016 | 0.0005 | 0.163 | 83.1% | 50.2% | 1.88% | 0.9937 | 0.42 pp |
| **Wall_Lizard** | 8 yr | 0.0150 | 0.0069 | 0.0038 | 0.0018 | 0.250 | 44.6% | 36.7% | 3.51% | 0.9941 | 0.30 pp |
| **Honeybee** | 8 yr | 0.0359 | 0.0162 | 0.0102 | 0.0066 | 0.283 | 53.8% | 67.3% | 7.28% | 0.9818 | 6.20 pp |
| **Green_Anole** | 7.2 yr | 0.0131 | 0.0072 | 0.0035 | 0.0017 | 0.265 | 50.4% | 36.2% | 2.83% | 0.9957 | 1.12 pp |
| **Fire_Ant** | 7 yr | 0.0160 | 0.0117 | 0.0074 | 0.0045 | 0.461 | 48.4% | 51.7% | 6.30% | 0.9259 | 0.02 pp |
| **House_Mouse** | 4 yr | 0.0176 | 0.0127 | 0.0064 | 0.0023 | 0.364 | 42.2% | 53.1% | 2.19% | 0.9973 | 0.64 pp |
| **Turquoise_Killifish** | 0.6 yr | 0.0270 | 0.0178 | 0.0105 | 0.0037 | 0.388 | 62.3% | 56.9% | 4.08% | 0.9854 | 1.62 pp |
| **Fruit_Fly** | 0.25 yr | 0.0149 | 0.0082 | 0.0040 | 0.0018 | 0.267 | 65.6% | 70.0% | 3.74% | 0.9909 | 1.42 pp |
| **Malaria_Mosquito** | 0.15 yr | 0.0227 | 0.0109 | 0.0035 | 0.0020 | 0.156 | 59.9% | 76.5% | 4.23% | 0.9920 | 0.18 pp |
| **Yellow_Fever_Mosquito** | 0.15 yr | 0.0068 | 0.0057 | 0.0023 | 0.0012 | 0.335 | 40.4% | 18.3% | 6.25% | 0.9951 | 0.12 pp |
| **Roundworm_Celegans** | 0.06 yr | 0.0132 | 0.0077 | 0.0045 | 0.0018 | 0.343 | 39.7% | 11.4% | 4.83% | 0.9006 | 4.11 pp |

## 2. Toric Evolutionary Slowdown Clock (Beyond-Markov-2 Drift Rate per 100 Myr)

| Species Pair | Clade | Divergence (Myr) | cos(Beyond-M2) | Toric Dist (1−cos) | Drift / 100 Myr |
|---|---|---:|---:|---:|---:|
| **Aldabra_Tortoise** × **Galapagos_Tortoise** | Testudines (Giant Tortoises) | 38 | 0.9972 | 0.0028 | **0.0075** |
| **Aldabra_Tortoise** × **Desert_Tortoise** | Testudines (Tortoises) | 55 | 0.9926 | 0.0074 | **0.0135** |
| **Aldabra_Tortoise** × **Painted_Turtle** | Testudines (Cryptodira) | 82 | 0.9805 | 0.0195 | **0.0238** |
| **Aldabra_Tortoise** × **American_Alligator** | Archelosauria | 255 | 0.8901 | 0.1099 | **0.0431** |
| **Ocean_Quahog** × **Pacific_Oyster** | Bivalvia (507y vs 20y) | 460 | 0.7008 | 0.2992 | **0.0650** |
| **Honeybee** × **Black_Garden_Ant** | Aculeata (Bee vs Ant) | 160 | 0.8776 | 0.1224 | **0.0765** |
| **Black_Garden_Ant** × **Fire_Ant** | Formicidae (Lasius vs Fire Ant) | 100 | 0.9031 | 0.0969 | **0.0969** |
| **Aldabra_Tortoise** × **Softshell_Turtle** | Testudines (Deep) | 160 | 0.8112 | 0.1888 | **0.1180** |
| **Wall_Lizard** × **Green_Anole** | Squamata (Lizards) | 165 | 0.8034 | 0.1966 | **0.1192** |
| **Rougheye_Rockfish** × **Turquoise_Killifish** | Teleostei (205y vs 0.6y) | 115 | 0.8591 | 0.1409 | **0.1226** |
| **Fruit_Fly** × **Malaria_Mosquito** | Diptera (Fly vs Anopheles) | 260 | 0.6225 | 0.3775 | **0.1452** |
| **Naked_Mole_Rat** × **House_Mouse** | Rodentia (Mole-Rat vs Mouse) | 70 | 0.8510 | 0.1490 | **0.2129** |
| **Human_T2T** × **Blue_Whale** | Boreoeutheria (Human-Whale) | 94 | 0.6935 | 0.3065 | **0.3261** |
| **Human_T2T** × **House_Mouse** | Euarchontoglires | 87 | 0.7113 | 0.2887 | **0.3318** |
| **Human_T2T** × **Greater_Mouse_Eared_Bat** | Boreoeutheria (Human-Bat) | 94 | 0.6626 | 0.3374 | **0.3589** |
| **Human_T2T** × **African_Elephant** | Placentalia (Human-Elephant) | 99 | 0.5546 | 0.4454 | **0.4499** |
| **Malaria_Mosquito** × **Yellow_Fever_Mosquito** | Culicidae (Anopheles vs Aedes) | 150 | 0.2679 | 0.7321 | **0.4881** |

## 3. Theorems L5d, L8b & L8c — Metazoan Structural Invariants (Zero-RNG Exact DP)

- **THM-L8b (447/512 Isotropy Relaxation):** Exact ceiling `87.30%`; 29 non-fly species relax at `d=30 bp` to `48.3%–85.9%`, while `Fruit_Fly` dips to `29.9%` at `d=10 bp`.
- **THM-L8c (Mammalian Palindromic Depletion):** `n1=6` mammals (med `2.05%`) vs `n2=24` non-mammals (med `3.40%`), `U1=16/144`, `r_rb=-0.778`, exact DP `p=2.14e-03`.
- **THM-L5d (Amniote Hurst Bifurcation):** `n1=13` full-chr amniotes (>100 Mb, med `0.850`) vs `n2=14` non-amniotes (med `0.602`), `r_rb=+0.857`, exact DP `p=3.72e-05`.
- **THM-L8d (Sauropsid Markov-2 Fidelity Law):** `n=8` sauropsids, `rho_S(cos_M2, Life)=+1.000`, exact permutation `p=4.96e-05`.
- **THM-L8e (Mammalian Markov-2 Overshoot Inversion):** `n1=6` mammals (med `+5.81` pp) vs `n2=11` non-mammal vertebrates (med `-9.70` pp), `r_rb=+1.000`, exact DP `p=1.62e-04`.
- **THM-L5e (COSMIC Clade Bifurcations):** Amniotes vs non-amniotes `SBS7a` UV `r_rb=+0.920` (exact DP `p=1.33e-06`), `SBS13` chiral share `r_rb=-0.562` (DP `p=7.92e-03`); Mammals vs Sauropsids (`SBS2` & `TCW` O/E) `r_rb=+1.000` (exact DP `p=6.66e-04`).

