# Final-Paper Checkpoint Manifest

Status: the exact final-paper checkpoint files are not present locally. Do not upload the available July checkpoints as paper models.

All reported runs use split SHA-256:

    dc26fae17e54fd2ad41a9e10353b3da3e0aacf3144b64f6ee62e8341b4360555

## Pfam-Augmented SetNet v5

| Seed | Validation-selected SetNet weight | Expected checkpoint SHA-256 |
| ---: | ---: | --- |
| 20260810 | 0.1 | bacada398caf49604b8b0c37842775ec2d517c7a0f244fa849e313dca3e68729 |
| 20260811 | 0.1 | 55adc3a7aecb9d1f65c916c2c31c9230358cb5dd65865f0566cad19595d67201 |
| 20260812 | 0.1 | da9b8ec8fc47d53b9abdc56eedbfd2468ba87ee3fe69980ad156d64c6d4016e6 |
| 20260813 | 0.2 | c537fd32ad380e0a8b7c9531530fcf7d9a3f73366a904c545666bd41462647dd |
| 20260814 | 0.4 | 9320902be03fb452ab907fa5b5ae6334ae632f72cb941715b85e7e43dbcc7153 |

## Weighted Pfam v6

| Seed | Expected checkpoint SHA-256 |
| ---: | --- |
| 20260810 | 75abede68918a7170ffd6b02881ce4bbdacc36ec028ee046c50a7ee1ae4fc899 |
| 20260811 | b05ab72dc19a1b9f06c6715eecaf9d0262b895ef4c63a32101beec15d0174677 |
| 20260812 | 215dd63ef93982389d826558497b95bbb3db157cfa24844541a9453523729b37 |
| 20260813 | 6b765c5fe7ef210f25dbaeeca07a7cff4b14afe8b6fdcf0c84526eb073af34e5 |
| 20260814 | 5e3265698ddc99c51ba263e5c55a008a962f288258fe0d6dd647bbbfd2be8553 |

## Required Before Hugging Face Publication

Recover all ten checkpoints from the original DGX archive and verify every hash above. Also recover the exact v5/v6 training and inference source, phase-1 provenance, model configuration, 2,218-token training-only Pfam vocabulary including PAD/UNK mapping, preprocessing specification, frozen split manifest, and artifact licensing.

All five seeds must be released. Selecting a seed using held-out test performance would be post-hoc model selection.
