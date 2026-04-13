import os

new_bib = """

@article{arnold_are_2000,
  title={Are tropical fungal endophytes hyperdiverse?},
  author={Arnold, A Elizabeth and Maynard, Zanda and Gilbert, Gregory S and Coley, Phyllis D and Kursar, Thomas A},
  journal={Ecology letters},
  volume={3},
  number={4},
  pages={267--274},
  year={2000},
  publisher={Wiley Online Library}
}

@book{zotz_plants_2016,
  title={Plants on plants-the biology of vascular epiphytes},
  author={Zotz, Gerhard},
  year={2016},
  publisher={Springer}
}

@article{vega_fungal_2008,
  title={Fungal endophytes in coffee plants from Colombia, Hawai'i, and Puerto Rico},
  author={Vega, Fernando E and Posada, Flavia and Aime, M Catherine and Pava-Ripoll, M and Infante, Francisco and Rehner, Stephen A},
  journal={Fungal Diversity},
  volume={31},
  pages={122},
  year={2008}
}

@article{dearnaley_orchid_2012,
  title={Orchid mycorrhizas: molecular ecology, physiology, evolution and conservation aspects},
  author={Dearnaley, John DW and Martos, Florent and Selosse, Marc-Andr{\\'e}},
  journal={Fungal associations},
  pages={207--230},
  year={2012},
  publisher={Springer}
}
"""

with open('/Users/mathys/Documents/Etudes/Stage_projet_LOT/CRBE/Rapport/Academique/BIBLIO_Stage_LOT_CRBE.bib', 'a') as f:
    f.write(new_bib)
