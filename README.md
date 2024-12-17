# MaxDGram-MapPigments

Map pigments in hyperspectral images (HSI) using endmembers and spectral angle mapping (SAM). Endmembers are extracted directly from the image using the Maximum Distance (MaxD) and Gram Matrix methods, enabling precise pigment identification and analysis.

## Project Structure

```
MaxDGram-MapPigments/
├── src/
│   ├── MaxD_Gram.py  -- Extracts end members using MaxD and estimates the material diversity using Gram matrix.
│   ├── spectral_tools.py -- Simple spectral tools, starting with spectral angle calculation.
├── pigment_map_MaxD.ipynb -- The notebook for analyzing and visualizing pigment distribution using maximum distance.
├── pigment_map_MaxD.py -- Script for mapping pigments using the MaxDGram algorithm
├── MaxD_Gram_pigmentmapping_pipeline.md -- The pipeline of the MaxDGram-MapPigments project, detailing the steps and processes involved.
└── README.md
```

## Dependencies
```bash
pip install numpy matplolib spectral
```


## Citation

If you find this repository useful in your research, please consider the following citation.

```bib
@article{kleynhans2020towards,
  title={Towards automatic classification of diffuse reflectance image cubes from paintings collected with hyperspectral cameras},
  author={Kleynhans, Tania and Messinger, David W and Delaney, John K},
  journal={Microchemical Journal},
  volume={157},
  pages={104934},
  year={2020},
  publisher={Elsevier}
}
```