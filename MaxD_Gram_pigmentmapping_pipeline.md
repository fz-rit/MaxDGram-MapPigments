# MaxD Pipeline

```mermaid
flowchart TD
    Start[Start: Pigment Mapping with MaxDist-Gram and SAM]
    Start --> A["Select a relatively large number of endmembers (e.g., 15)"]
    A --> B[Calculate vector lengths of each pixel]
    B --> C[Select pixels with largest and smallest distances as first two endmembers]
    C --> D[Project data onto subspaces orthogonal to difference vectors]
    D --> E[Iteratively find subsequent endmembers]
    E --> F[Convert endmembers into a Gram Matrix]
    F --> G[Calculate determinant to identify linearly independent sets]
    G --> H["Plot volumes (square roots of determinants) of these sets"]
    H --> I[Determine number of distinct endmember spectra from plot]
    I --> J["Create classification maps using Spectral Angle Mapper (SAM)"]
    J --> K[Assign each pixel to endmember with smallest spectral angle]
    K --> L{Are histograms of angles bi-modal?}
    L -->|Yes| M[Divide endmember into sub-classes]
    M --> N[Create new endmembers by averaging spectra in secondary mode]
    N --> End[End]
    L -->|No| End
