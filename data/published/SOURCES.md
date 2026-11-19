# Source documents

Every external number used in this study is listed, with its source and location, in
`benchmarks.csv`, `epc_measured_boards.csv` and `jlcpcb_stackups.csv`. The source
documents themselves are copyrighted and are **not** redistributed here. To re-run the
digitisation step (`scripts/digitise_epc.py`), obtain them from the publishers and place
them in `data/published/sources/` under these file names:

| File name expected | Document |
|---|---|
| `Reusch_APEC2013_slides.pdf` | D. Reusch and J. Strydom (EPC), "Understanding the effect of PCB layout on circuit performance in a high frequency gallium nitride based point of load converter", APEC 2013 presentation slides |
| `EPC_WP010_Optimizing_PCB_Layout.pdf` | Efficient Power Conversion, White Paper WP010, "Optimizing PCB Layout" (2019), epc-co.com |
| `EPC2015C_datasheet.pdf` | Efficient Power Conversion, EPC2015C datasheet (rev. April 2021), epc-co.com |
| `EPC_HighPerf_Layout_webinar.pdf` | Efficient Power Conversion, "High-Performance Layout Techniques to Maximize GaN Device Performance", How2GaN webinar slides (2021) |
| `Dubois_TPEL2018_early.pdf` | A. Letellier, M. Dubois, J. P. F. Trovao, H. Maher, IEEE Trans. Power Electron., DOI 10.1109/TPEL.2018.2826920 |

Also used (not needed to run the code): D. Reusch and J. Strydom, IEEE Trans. Power
Electron. 29(4), 2014, DOI 10.1109/TPEL.2013.2266103; JLCPCB stack-up and capability pages
(https://jlcpcb.com/impedance, https://jlcpcb.com/capabilities/pcb-capabilities).

The digitised values (`digitised_apec13.csv`, `epc_fig7b_digitised.csv`) are included, so
all downstream results can be reproduced without the source documents.
